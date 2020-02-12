#pragma once
#include <fstream>
#include "../../Maths/Spaces.h"

template <typename Space>
void WriteVelocityHandle(std::fstream& file);

template <>
void WriteVelocityHandle<Space2>(std::fstream& file)
{
  file << "Vx;Vy;";
}

template <>
void WriteVelocityHandle<Space3>(std::fstream& file)
{
  file << "Vx;Vy;Vz;";
}

template <typename Space>
void ReadSample(FILE* file, typename Space::Scalar& time, typename Space::Vector& v, typename Space::Tensor& sigma);

template <>
void ReadSample<Space2>(FILE* file, Space2::Scalar& time, Space2::Vector& v, Space2::Tensor& sigma)
{
  int itemsCount = fscanf(file, "%lf;%lf;%lf;%lf;%lf;%lf;\n",
    &time, &(sigma.xx), &(sigma.yy), &(sigma.xy), &(v.x), &(v.y));
  if (itemsCount != 6) std::cerr << "Can`t read line\n";
}

template <>
void ReadSample<Space3>(FILE* file, Space3::Scalar& time, Space3::Vector& v, Space3::Tensor& sigma)
{
  int itemsCount = fscanf(file, "%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;\n",
    &time, &(sigma.xx), &(sigma.yy), &(sigma.zz), &(sigma.xy), &(sigma.yz), &(sigma.xz),
    &(v.x), &(v.y), &(v.z));
  if (itemsCount != 9) std::cerr << "Can`t read line\n";
}

template <typename Space>
void WriteSample(std::fstream& velocityFile, std::fstream& pressureFile,
  typename Space::Scalar time, const typename Space::Vector& v, const typename Space::Tensor& sigma);

template <>
void WriteSample<Space2>(std::fstream& velocityFile, std::fstream& pressureFile,
  Space2::Scalar time, const Space2::Vector& v, const Space2::Tensor& sigma)
{
  if(velocityFile.is_open())
    velocityFile << v.x << ";" << v.y << ";";
  if(pressureFile.is_open())
    pressureFile << -(sigma.xx + sigma.yy) / Space2::Scalar(2.0) << ";";
}

template <>
void WriteSample<Space3>(std::fstream& velocityFile, std::fstream& pressureFile,
  Space3::Scalar time, const Space3::Vector& v, const Space3::Tensor& sigma)
{
  if(velocityFile.is_open())
    velocityFile << v.x << ";" << v.y << ";" << v.z << ";";
  if(pressureFile.is_open())
    pressureFile << -(sigma.xx + sigma.yy + sigma.zz) / Space3::Scalar(3.0) << ";";
}
