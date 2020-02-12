#pragma once

#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

bool ReplaceSubstring(std::string& str, const std::string& from, const std::string& to)
{
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
      return false;
  str.replace(start_pos, from.length(), to);
  return true;
}

template <typename T>
size_t RemoveDuplicates(std::vector<T>* const v)
{
  std::sort(v->begin(), v->end());
  typename std::vector<T>::iterator it;
  it = std::unique(v->begin(), v->end()); 
  v->resize(std::distance(v->begin(),it));
  return v->size();
}

std::vector<std::string>& Split(const std::string &s, char delim, std::vector<std::string> &elems)
{
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
  {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> Split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  Split(s, delim, elems);
  return elems;
}

template<typename T>
std::set<T> StringToSetOfInt(const std::string& s, char delim = ' ')
{
  std::set<T> res;

  std::vector<std::string> v = Split(s, delim);
  for (size_t i = 0; i < v.size(); ++i)
  {
    T x = static_cast<T>(std::stoi(v[i]));
    res.insert(x);
  }

  return res;
}

bool EndWith(const std::string& fullString, const std::string& ending)
{
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  }
  else {
    return false;
  }
}

std::string AddExtensionToFileName(const std::string& fileName, const std::string& extension)
{
  std::string result = fileName;

  if (!EndWith(result, extension))
   result += extension;

  return result;
}

template <typename T>
void swap_endian(T& pX)
{
  // should static assert that T is a POD

  char& raw = reinterpret_cast<char&>(pX);
  std::reverse(&raw, &raw + sizeof(T));
}

template <typename T>
T swap_endian_copy(T pX)
{
  swap_endian(pX);
  return pX;
}

namespace IO
{
  enum FileType 
  {
    Ascii = 0, Binary = 1
  };

  template <typename T>
  void Write(std::fstream& file, T x)
  {
    file.write((const char *)(&x), sizeof(T));
  }

  template <typename T>
  void Write(std::fstream& file, T* const x, size_t elementsCount)
  {
    file.write((const char *)x, sizeof(T) * elementsCount);
  }

  template <typename T>
  void Write(std::fstream& file, T x, T y)
  {
    Write<T>(file, x);
    Write<T>(file, y);
  }

  template <typename T>
  void Write(std::fstream& file, T x, T y, T z)
  {
    Write<T>(file, x);
    Write<T>(file, y);
    Write<T>(file, z);
  }

  template <typename T>
  void WriteVector(std::fstream& file, const std::vector<T>& v)
  {
    Write(file, v.data(), v.size());
  }

  template <typename T>
  void Read(std::fstream& file, T& x)
  {
    file.read((char*)&x, sizeof(T));
  }

  template <typename T>
  void Read(std::fstream& file, T* x, size_t elementsCount)
  {
    if (x && elementsCount)
      file.read((char*)x, sizeof(T) * elementsCount);
  }

  template <typename T>
  void Read(std::fstream& file, T& x, T& y)
  {
    Read<T>(file, x);
    Read<T>(file, y);
  }

  template <typename T>
  void Read(std::fstream& file, T& x, T& y, T& z)
  {
    Read<T>(file, x);
    Read<T>(file, y);
    Read<T>(file, z);
  }

  template <typename T>
  void Save(std::fstream& file, FileType fileType, const std::vector<T>& v)
  {
    switch (fileType)
    {
      case Ascii:
        file << v.size() << std::endl;
        for (size_t index = 0; index < v.size(); ++index)
          file << v[index] << " ";
        if (!v.empty()) 
          file << std::endl;
      break;
      case Binary:
        Write(file, v.size());
        WriteVector(file, v);
      break;
    }
  }

  template <typename T>
  void Load(std::fstream& file, FileType fileType, std::vector<T>* const v)
  {
    size_t size;
    switch (fileType)
    {
      case Ascii:
        file >> size;
        v->resize(size);
        for (size_t index = 0; index < v->size(); ++index)
          file >> (*v)[index];
      break;
      case Binary:
        Read(file, size);
        v->resize(size);
        Read(file, v->data(), size);
      break;
    }
  }
}
