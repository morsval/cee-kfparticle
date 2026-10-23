#ifndef CEE_KFP_INPUT_FILES_H
#define CEE_KFP_INPUT_FILES_H
#include "TList.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace CeeV0Analysis {
namespace InputFiles {
inline std::string Trim(const std::string &input) {
  const std::string whitespace = " \t\r\n";
  const std::size_t first = input.find_first_not_of(whitespace);
  if (first == std::string::npos)
    return std::string();
  const std::size_t last = input.find_last_not_of(whitespace);
  return input.substr(first, last - first + 1u);
}

inline bool EndsWith(const std::string &value, const std::string &suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) ==
             0;
}

inline std::vector<std::string> ResolveInputFiles(const std::string &input,
                                                  int maximumFiles) {
  std::vector<std::string> files;
  if (EndsWith(input, ".root")) {
    files.push_back(input);
  } else if (EndsWith(input, ".list") || EndsWith(input, ".txt")) {
    std::ifstream stream(input.c_str());
    if (!stream)
      throw std::runtime_error("cannot read input list: " + input);
    std::string line;
    while (std::getline(stream, line)) {
      line = Trim(line);
      if (line.empty() || line[0] == '#')
        continue;
      if (!EndsWith(line, ".root"))
        throw std::runtime_error("invalid ROOT path in list: " + line);
      files.push_back(line);
    }
  } else {
    TSystemDirectory directory("ceeV0Input", input.c_str());
    TList *entries = directory.GetListOfFiles();
    if (entries) {
      TIter next(entries);
      while (TSystemFile *entry = static_cast<TSystemFile *>(next())) {
        const std::string name = entry->GetName();
        if (!entry->IsDirectory() && EndsWith(name, ".root"))
          files.push_back(input + "/" + name);
      }
    }
    std::sort(files.begin(), files.end());
  }
  if (maximumFiles > 0 && files.size() > static_cast<std::size_t>(maximumFiles))
    files.resize(static_cast<std::size_t>(maximumFiles));
  return files;
}

} // namespace InputFiles
} // namespace CeeV0Analysis
#endif
