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

int main(int argc, char** argv)
{
  if (argc < 2)
    return std::cerr << "Error: first argument must be random seed\n", EXIT_FAILURE;

  if (argc < 3)
    return std::cerr << "Error: second arguiment must be recom target\n", EXIT_FAILURE;

  std::int64_t seed = std::atoll(argv[1]);
  float target_records_per_recom = std::atof(argv[2]);

  savvy::reader in(argc > 3 ? argv[3] : "/dev/stdin");
  if (!in)
    return std::cerr << "Error: could not open input file\n", EXIT_FAILURE;

  if (in.samples().empty())
   return std::cerr << "Error: no samples in input file\n", EXIT_FAILURE;

  //std::size_t n_haps = in.samples().size() * 2; // Assuming diploid



  auto ids = in.samples();
  for (std::size_t i = 0; i < ids.size(); ++i)
    ids[i] = std::to_string(i);

  savvy::writer out("/dev/stdout", savvy::file::format::bcf, in.headers(), ids, 0);

  double recom_prob = 0.5 / target_records_per_recom; // Using 0.5 because two haps are switched per event.
  std::cerr << "Recom prob: " << recom_prob << std::endl;
  //std::default_random_engine prng(seed);
  std::mt19937_64 prng(seed);
  std::geometric_distribution<std::int64_t> geom_dist(recom_prob);

  std::vector<std::int8_t> gt, gt_shuffled;
  savvy::variant rec;
  if (!(in >> rec))
    return std::cerr << "Error: empty VCF\n", EXIT_FAILURE;
  rec.get_format("GT", gt);

  std::vector<std::size_t> random_hap_idx;
  random_hap_idx.reserve(gt.size());
  for (std::size_t i = 0; i < gt.size(); ++i)
  {
    if (!savvy::typed_value::is_end_of_vector(gt[i]))
      random_hap_idx.emplace_back(i);
  }


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

  do //while (in >> rec)
  {
    std::int64_t remaining_hap_positions_in_row = random_hap_idx.size();

    while (true)
    {
      if (geom_draw >= remaining_hap_positions_in_row)
      {
        geom_draw -= remaining_hap_positions_in_row;
        break;
      }
      else
      {
        std::size_t hap1 = *random_hap_idx_it;
        increment_rand(random_hap_idx_it, random_hap_idx, prng);

        std::size_t hap2 = *random_hap_idx_it;
        increment_rand(random_hap_idx_it, random_hap_idx, prng);

        // if (hap1 != hap2) This will never happen with shuffle approach
        {
          ++hap_switch_cnts[hap1];
          ++hap_switch_cnts[hap2];
          std::swap(hap_mapping[hap1], hap_mapping[hap2]);
        }
      
        remaining_hap_positions_in_row -= 1 + geom_draw;
        geom_draw = geom_dist(prng);
        if (geom_draw < 0ll) throw std::runtime_error("should not be negative");
      }
    }

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
