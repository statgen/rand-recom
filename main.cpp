/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <cstdlib>
#include <random>
#include <algorithm>

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include "prog_args.hpp"

template<typename Prng>
void increment_rand(std::vector<std::size_t>::iterator& it, std::vector<std::size_t>& vec, Prng& r)
{
  ++it;
  if (it == vec.end())
  {
    std::shuffle(vec.begin(), vec.end(), r);
    it = vec.begin();
  }
}

void defer_rand(std::vector<std::size_t>::iterator& it, std::vector<std::size_t>& vec)
{
  if (it + 1 == vec.end())
    std::swap(*it, *(vec.begin()));
  else
    std::swap(*it, *(it + 1));
}

int main(int argc, char** argv)
{
  prog_args args;
  if (!args.parse(argc, argv))
  {
    args.print_usage(std::cerr);
    return EXIT_FAILURE;
  }

  if (args.help_is_set())
  {
    args.print_usage(std::cout);
    return EXIT_SUCCESS;
  }

  if (args.version_is_set())
  {
    std::cout << "rand-recom v" << RAND_RECOM_VERSION << std::endl;
    return EXIT_SUCCESS;
  }

  savvy::reader in(args.input_path());
  if (!in)
    return std::cerr << "Error: could not open input file\n", EXIT_FAILURE;

  if (in.samples().empty())
   return std::cerr << "Error: no samples in input file\n", EXIT_FAILURE;

  auto ids = in.samples();
  for (std::size_t i = 0; i < ids.size(); ++i)
    ids[i] = std::to_string(i);

  savvy::writer out(args.output_path(), args.output_format(), in.headers(), ids, args.output_compression_level());

  std::vector<std::int8_t> gt, gt_shuffled;
  savvy::variant rec;
  if (!(in >> rec))
    return std::cerr << "Error: empty VCF\n", EXIT_FAILURE;
  rec.get_format("GT", gt);

  double recom_prob = 0.5 / args.target_segment_length(); // Using 0.5 because two haps are switched per event.
  std::cerr << "Recom prob: " << recom_prob << std::endl;

  std::mt19937_64 prng(args.seed());
  std::geometric_distribution<std::int64_t> geom_dist(recom_prob);

  std::vector<std::size_t> random_hap_idx;
  random_hap_idx.reserve(gt.size());
  for (std::size_t i = 0; i < gt.size(); ++i)
  {
    if (!savvy::typed_value::is_end_of_vector(gt[i]))
      random_hap_idx.emplace_back(i);
  }

  std::vector<std::size_t> non_eov_mapping = random_hap_idx;
  std::vector<std::size_t> hap_mapping(gt.size());
  std::iota(hap_mapping.begin(), hap_mapping.end(), 0);
  std::vector<std::size_t> hap_switch_cnts(gt.size());

  std::shuffle(random_hap_idx.begin(), random_hap_idx.end(), prng);
  auto random_hap_idx_it = random_hap_idx.begin();

  for (std::size_t i = 0; i < gt.size(); ++i)
  {
    if (!savvy::typed_value::is_end_of_vector(gt[i]))
      std::swap(hap_mapping[i], hap_mapping[*(random_hap_idx_it++)]);
  }

  std::shuffle(random_hap_idx.begin(), random_hap_idx.end(), prng);
  random_hap_idx_it = random_hap_idx.begin();

  std::int64_t geom_draw = geom_dist(prng);
  std::size_t bp_pos = geom_draw / random_hap_idx.size();
  std::size_t hap_pos = geom_draw % random_hap_idx.size();

  do
  {
    while (bp_pos < rec.pos())
    {
      std::size_t h1;
      if (args.uniform())
      {
        h1 = *random_hap_idx_it;
        increment_rand(random_hap_idx_it, random_hap_idx, prng);
      }
      else
      {
        h1 = non_eov_mapping[hap_pos];
        if (*random_hap_idx_it == h1)
          increment_rand(random_hap_idx_it, random_hap_idx, prng); //defer_rand(random_hap_idx_it, random_hap_idx);
      }

      std::size_t h2 = *random_hap_idx_it;
      increment_rand(random_hap_idx_it, random_hap_idx, prng);

      ++hap_switch_cnts[h1];
      ++hap_switch_cnts[h2];
      std::swap(hap_mapping[h1], hap_mapping[h2]);

      geom_draw = geom_dist(prng) + 1; // +1 increments hap_pos

      if (hap_pos + geom_draw < random_hap_idx.size())
      {
        hap_pos += geom_draw;
      }
      else
      {
        assert(hap_pos <= random_hap_idx.size());
        assert(geom_draw >= random_hap_idx.size() - hap_pos);
        geom_draw -= random_hap_idx.size() - hap_pos;
        bp_pos += geom_draw / random_hap_idx.size();
        hap_pos = geom_draw % random_hap_idx.size();
      }
    }
#if 0
    // This is for testing chunked panel creation.
    if (rec.pos() >= 3400001)
    {
      for (auto it = hap_mapping.begin(); it != hap_mapping.end(); ++it)
        std::cout << *it << "\n";
      return EXIT_FAILURE;
    }
#endif
    rec.get_format("GT", gt);
    if (gt.size() != hap_mapping.size())
      return std::cerr << "Error: inconsistent ploidy\n", EXIT_FAILURE;

    gt_shuffled.resize(gt.size());
    for (std::size_t i = 0; i < gt.size(); ++i)
      gt_shuffled[hap_mapping[i]] = gt[i];

    rec.set_format("GT", gt_shuffled);

    out << rec;
  } while (in >> rec);

  if (in.bad() || !out.good())
    return std::cerr << "Error: I/O failure\n", EXIT_FAILURE;


  for (std::size_t i = 0; i < gt.size(); ++i)
  {
    if (savvy::typed_value::is_end_of_vector(gt[i]) && hap_switch_cnts[i] != 0)
      return std::cerr << "Error: switch occurred at EOV position in last variant (make sure all variants have consistent ploidy)\n", EXIT_FAILURE;
  }

  for (std::size_t i = 0; i < in.samples().size(); ++i)
    std::cerr << in.samples()[i] << "\t" << hap_switch_cnts[i*2] << "\t" << hap_switch_cnts[i*2+1] << "\n";
  return EXIT_SUCCESS;
}
