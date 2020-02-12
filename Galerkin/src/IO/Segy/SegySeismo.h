#pragma once

#include <string>
#include <vector>

struct segy_bin_header_data
{
  typedef unsigned int uint32;
  typedef unsigned short uint16;

  uint32 job_id;
  uint32 line_num;
  uint32 reel_num;
  uint16 num_of_traces_per_record;
  uint16 num_of_auxiliary_traces_per_record;
  uint16 sample_interval_reel;
  uint16 sample_interval;
  uint16 samples_per_trace_reel;
  uint16 samples_per_trace;
  // 	Data sample format code: 1 = IBM floating point (4 bytes) 2 = fixed point (4 bytes)
  //  3 = fixed point (2 bytes) 4 = fixed point with gain code (4 bytes), 5 - IEEE floating point.
  uint16 data_sample_format;
  // Skip other data
  uint16 other[17];
  // Reserve
  uint16 reserve[170];
};

struct segy_trace_header
{
  typedef unsigned int uint32;
  typedef unsigned short uint16;

  uint32 trace_seq_num_line;
  uint32 trace_seq_num_reel;
  uint32 field_record_num;
  uint32 trace_num_reel;

  uint16 other1[6];

  uint16 trace_id_code;

  uint16 num_of_verticaly_summed_traces;
  uint16 num_of_horizotally_summed_traces;
  uint16 data_use;
  uint32 distance_from_source;

  uint16 other2[15];

  uint32 source_x;
  uint32 source_y;
  uint32 receiver_x;
  uint32 receiver_y;
  // Coordinate units: 1 = length (meters or feet) 2 = seconds of arc
  uint16 units_id;

  uint16 other3[12];
  uint16 num_of_samples;
  uint16 sample_interval; // in microseconds

  uint16 other4[31];
  uint16 reserve[30];
};


template <typename Scalar>
class Seismogramm
{
public:
  typedef Scalar Sample;
  typedef std::vector<Sample> Trace;
  typedef size_t IndexType;
  typedef unsigned int uint32;
  typedef unsigned short uint16;

  Seismogramm() {}

  void LoadSegY(const std::string& path, std::vector<Scalar>& times);
  void SaveSegY(const std::string& path, const std::vector<Scalar>& times);
  void AddValue(const Sample& value, IndexType detectorIndex);

  std::vector<Trace> data;
  struct segy_bin_header_data header_data;

  std::vector<struct segy_trace_header> trace_header_data;

private:
  std::vector<Scalar> interpolationIntervals;
  void swap_header_endian();
  void swap_trace_header_endian(struct segy_trace_header * ptr_header);
  void swap_all_trace_headers_endian();
  void swap_data_endian();
};

enum SeismoType
{
  SEG_Y, CSV
};

// ////////////////////////////////////////////////////////
// Value getters:

template <typename Elastic, int dims>
struct ValueGetter
{
  typedef typename Elastic::ScalarType Scalar;
  virtual Scalar GetValue(const Elastic& elastic) const = 0;
};

// vx getter
template <typename Elastic, int dims>
struct VxGetter;

template <typename Elastic>
struct VxGetter<Elastic, 2>: public ValueGetter<Elastic, 2>
{
  typedef typename Elastic::ScalarType Scalar;
  Scalar GetValue(const Elastic& elastic) const
  {
    return elastic.values[0];
  }
};

template <typename Elastic>
struct VxGetter<Elastic, 3>: public ValueGetter<Elastic, 3>
{
  typedef typename Elastic::ScalarType Scalar;
  Scalar GetValue(const Elastic& elastic) const
  {
    return elastic.values[0];
  }
};

// vy getter
template <typename Elastic, int dims>
struct VyGetter;

template <typename Elastic>
struct VyGetter<Elastic, 2>: public ValueGetter<Elastic, 2>
{
  typedef typename Elastic::ScalarType Scalar;
  Scalar GetValue(const Elastic& elastic) const
  {
    return elastic.values[1];
  }
};

template <typename Elastic>
struct VyGetter<Elastic, 3>: public ValueGetter<Elastic, 3>
{
  typedef typename Elastic::ScalarType Scalar;
  Scalar GetValue(const Elastic& elastic) const
  {
    return elastic.values[1];
  }
};

// pressure getter
template <typename Elastic, int dims>
struct PressureGetter;

template <typename Elastic>
struct PressureGetter<Elastic, 2>: public ValueGetter<Elastic, 2>
{
  typedef typename Elastic::ScalarType Scalar;
  Scalar GetValue(const Elastic& elastic) const
  {
    return -(elastic.values[2] + elastic.values[3]) / Scalar(2.0);
  }
};

template <typename Elastic>
struct PressureGetter<Elastic, 3>: public ValueGetter<Elastic, 3>
{
  typedef typename Elastic::ScalarType Scalar;
  Scalar GetValue(const Elastic& elastic) const
  {
    return -(elastic.values[3] + elastic.values[4] + elastic.values[5]) / Scalar(3.0);
  }
};

// ////////////////////////////////////////////////////////////////////

template <typename Scalar, int dims>
class CombinedSeismogramm
{
public:
  std::vector<Scalar> times;
  std::vector<Seismogramm<Scalar> > seismogramms;
  Scalar interpolation_multiplier;
  typedef size_t IndexType;
  typedef unsigned int uint32;
  typedef unsigned short uint16;

  struct Elastic
  {
    typedef Scalar ScalarType;
    Scalar values[(dims + 1) * (dims + 2) / 2 - 1]; //WHAT THE FUCK
  };

  struct ComponentInfo
  {
    ComponentInfo(const std::string& path, ValueGetter<Elastic, dims>* getter):
    path(path), getter(getter)
    {
    }
    std::string path;
    ValueGetter<Elastic, dims>* getter;
  };

  std::vector<ComponentInfo> componentInfos;

  CombinedSeismogramm(Scalar interpolation_multiplier = 1) : interpolation_multiplier(interpolation_multiplier) {}

  void Load(SeismoType type, std::vector<std::string> paths, const CombinedSeismogramm<Scalar, dims> *compatibleSeismogram = 0);

  void Save(SeismoType type, std::vector<std::string> paths);
  void Save(SeismoType type);

  void AddValue(Scalar time, const Elastic& elastic, IndexType detectorIndex);

  void AddComponent(const std::string& path, ValueGetter<Elastic, dims>* getter);

  CombinedSeismogramm<Scalar, dims> &operator -=(const CombinedSeismogramm<Scalar, dims> & other);

private:
  std::vector<Scalar> interpolationIntervals;
  void interpolate_data_on_equal_time_intervals(Scalar time_interval);
};

#include "SegySeismo.inl"