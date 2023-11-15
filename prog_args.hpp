/*
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "getopt_wrapper.hpp"

class prog_args : public getopt_wrapper
{
private:
  std::string input_path_ = "/dev/stdin";
  std::string recom_map_path_;
  std::string output_path_ = "/dev/stdout";
  savvy::file::format output_format_ = savvy::file::format::bcf;
  int output_compression_level_ = 6;
  std::int64_t target_segment_length_ = 0;
  std::int64_t seed_ = 0;
  bool uniform_ = false;
  bool version_ = false;
  bool help_ = false;
public:
  prog_args() :
    getopt_wrapper("Usage: rand-recom [opts ...] input.{sav,bcf,vcf.gz} --seed <integer> --target-length <bp_length>",
      {
        {"help", "", 'h', "Print usage"},
        {"version", "", 'v', "Print version"},
        {"output", "<file>", 'o', "Output path (default: /dev/stdout)"},
        {"output-format", "<string>", 'O', "Output file format (bcf, sav, vcf.gz, ubcf, usav, or vcf; default: bcf)"},
        {"seed", "<integer>", 's', "Seed for pseudorandom number generator"},
        {"target-length", "<integer>", 't', "Target segment length in bp"},
        {"uniform", "", 'u', "Enforces uniform distribution of switch counts per haplotype"},
      })
  {
  }

  const std::string& input_path() const { return input_path_; }
  const std::string& output_path() const { return output_path_; }
  const std::string& recom_map_path() const { return recom_map_path_; }

  std::int64_t seed() const { return seed_; }
  std::int64_t target_segment_length() const { return target_segment_length_; }
  savvy::file::format output_format() const { return output_format_; }
  int output_compression_level() const { return output_compression_level_; }

  bool uniform() const { return uniform_; }
  bool version_is_set() const { return version_; }
  bool help_is_set() const { return help_; }

  bool parse(int argc, char** argv)
  {
    int long_index = 0;
    int opt = 0;
    while ((opt = getopt_long(argc, argv, short_opt_string_.c_str(), long_options_.data(), &long_index )) != -1)
    {
      char copt = char(opt & 0xFF);
      switch (copt)
      {
      case 'h':
        help_ = true;
        return true;
      case 'v':
        version_ = true;
        return true;
      case 'u':
        uniform_ = true;
        break;
      case 'o':
        output_path_ = optarg ? optarg : "";
        break;
      case 'O':
        {
          using fmt = savvy::file::format;
          std::string ot = optarg ? optarg : "";
          if (ot == "vcf")
          {
            output_format_ = fmt::vcf;
            output_compression_level_ = 0;
          }
          else if (ot == "vcf.gz")
          {
            output_format_ = fmt::vcf;
          }
          else if (ot == "bcf")
          {
            output_format_ = fmt::bcf;
          }
          else if (ot == "ubcf")
          {
            output_format_ = fmt::bcf;
            output_compression_level_ = 0;
          }
          else if (ot == "sav")
          {
            output_format_ = fmt::sav;
          }
          else if (ot == "usav")
          {
            output_format_ = fmt::sav;
            output_compression_level_ = 0;
          }
          else
          {
            std::cerr << "Invalid --output-format: " << ot << std::endl;
            return false;
          }
          break;
        }
      case 't':
        {
          target_segment_length_ = std::atoll(optarg ? optarg : "");
          if (target_segment_length_ <= 0)
            return std::cerr << "Error: invalid --target-length\n", false;
          break;
        }
      case 's':
        {
          seed_ = std::atoll(optarg ? optarg : "");
          break;
        }
      default:
        return false;
      }
    }

    int remaining_arg_count = argc - optind;

    if (remaining_arg_count == 1)
    {
      input_path_ = argv[optind];
    }
    else if (remaining_arg_count > 1)
    {
      std::cerr << "Too many arguments\n";
      return false;
    }

    if (!seed_)
      return std::cerr << "Error: must provide a non-zero --seed argument\n", false;

    //if (!target_segment_length_)
    //  return std::cerr << "Error: must provide a --target-length argument\n", false;

    return true;
  }
};
