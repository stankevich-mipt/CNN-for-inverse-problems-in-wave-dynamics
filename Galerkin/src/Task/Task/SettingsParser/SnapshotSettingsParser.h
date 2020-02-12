#pragma once

#include "ParserUtil.h"

template<typename Space>
struct SnapshotSettings
{
  SPACE_TYPEDEFS

  SnapshotSettings(): 
    timePeriod(std::numeric_limits<Scalar>::infinity()),
    framePeriod(IndexType(-1)),
    startTime(Scalar(0.0)),
    finishTime(std::numeric_limits<Scalar>::infinity()),
    startFrame(0),
    finishFrame(IndexType(-1))
  {}

  struct Data
  {
    Data():
      filename("snapshot<step>"),
      used(true),
      writeVelocity(true),
      writeTension(false),
      writeWaves(false)
    {
      Initialize();
    }

    void Initialize()
    {
      Initialize(Overload<Space>());
      origin = Vector::zeroVector();
      boxPoint1 = Vector::zeroVector();
      boxPoint2 = Vector::zeroVector();
    }

    std::string filename;
    Vector origin;
    Vector spacing;
    IndexVector resolution;
    Vector boxPoint1;
    Vector boxPoint2;
    bool used;
    bool writeVelocity;
    bool writeTension;
    bool writeWaves;

  private:
    void Initialize(Overload<Space2>)
    {
      spacing = Vector(-1.0, -1.0);
      resolution = IndexVector(-1, -1);
    }
    void Initialize(Overload<Space3>)
    {
      spacing = Vector(-1.0, -1.0, -1.0);
      resolution = IndexVector(-1, -1, -1);
    }
  } data;

  struct Mesh
  {
    Mesh():
      filename("out/<domain>_mesh<step>.vtk"),
      used(true)
    {}
    std::string filename;
    bool used;
  } mesh;

  struct Contacts
  {
    Contacts():
      filename("out/contacts[<domain>]<step>.vtk"),
      used(true)
    {
      const int MaxContactCount = 100;
      std::vector<IndexType> v(MaxContactCount);
      for (int i = 0; i < MaxContactCount; ++i)
      {
        v[i] = i;
      }
      faceTypesToDraw.insert(v.begin(), v.end());
    }
    std::string filename;
    bool used;
    std::set<IndexType> faceTypesToDraw;
  } contacts;

  struct Boundaries
  {
    Boundaries():
      fileName("out/boundaries[<domain>]<step>.vtk"),
      used(true)
    {
      const int MaxContactCount = 100;
      std::vector<IndexType> v(MaxContactCount);
      for (int i = 0; i < MaxContactCount; ++i)
      {
        v[i] = i;
      }
      faceTypesToDraw.insert(v.begin(), v.end());
    }
    std::string fileName;
    bool used;
    std::set<IndexType> faceTypesToDraw;
  } boundaries;

  struct CellInfos
  {
    CellInfos() :
      fileName("out/cellinfo[<domain>]<step>.vtk"),
      used(true)
    {}
    std::string fileName;
    bool used;
  } cellInfos;

  Scalar timePeriod;
  Scalar startTime;
  Scalar finishTime;

  IndexType framePeriod;
  IndexType startFrame;
  IndexType finishFrame;


  void Parse(TiXmlElement* snapshotElement);
};

template<typename Space>
struct DetectorsSettings
{
  DetectorsSettings():
    filename("out/<domain>_detectors<detector>.csv"),
    locationsFileName("meshes/detectorsLocations.txt"),
    used(true), samplesPeriod(1)
  {}
  std::string filename;
  std::string locationsFileName;
  bool used;
  typename Space::IndexType samplesPeriod;

  void Parse(TiXmlElement* settingsElement);
};

template<typename Space>
void SnapshotSettings<Space>::Parse(TiXmlElement *snapshotElement)
{
  // periods
  TiXmlElement *periodElement = snapshotElement->FirstChildElement("Period");
  if(periodElement)
  {
    ParseScalar(periodElement, "time",     &timePeriod);
    ParseUnsigned(periodElement, "frames", &framePeriod);
  }

  if (timePeriod == std::numeric_limits<Scalar>::infinity() && framePeriod == IndexType(-1))
  {
    std::cerr << "There is no period attribute";
    throw;
  }

  TiXmlElement *intervalElement = snapshotElement->FirstChildElement("Interval");
  if(intervalElement)
  {
    ParseScalar(intervalElement, "startTime", &startTime);
    ParseScalar(intervalElement, "finishTime", &finishTime);

    ParseUnsigned(intervalElement, "startFrame", &startFrame);
    ParseUnsigned(intervalElement, "finishFrame", &finishFrame);
  }

  TiXmlElement* dataElement = snapshotElement->FirstChildElement("Data");
  if (!dataElement)
  {
    data.used = false;
  } else
  {
    ParseString(dataElement, "fileName",      &data.filename);
    ParseBool  (dataElement, "writeVelocity", &data.writeVelocity);
    ParseBool  (dataElement, "writeTension",  &data.writeTension);
    ParseBool  (dataElement, "writeWaves",    &data.writeWaves);

    TiXmlElement* boxElement = dataElement->FirstChildElement("Box");
    if (boxElement)
    {
      ParseVector(boxElement, "boxPoint1", &data.boxPoint1);
      ParseVector(boxElement, "boxPoint2", &data.boxPoint2);
    }

    // frame
    TiXmlElement* frameElement = dataElement->FirstChildElement("Frame");
    if (frameElement)
    {
      Vector center = Vector::zeroVector();
      Vector size   = Vector::one();
      Scalar margin = 0;

      ParseVector(frameElement, "center", &center);
      ParseVector(frameElement, "size"  , &size);
      ParseScalar(frameElement, "margin", &margin);

      data.boxPoint1 = center - size / Scalar(2.0) - Vector::one() * margin;
      data.boxPoint2 = center + size / Scalar(2.0) + Vector::one() * margin;
    }

    // resolution
    TiXmlElement* resolutionElement = dataElement->FirstChildElement("Resolution");
    if (resolutionElement)
    {
      ParseVector(resolutionElement, "resolution",  &data.resolution);

      data.origin  = data.boxPoint1;
      data.spacing = (data.boxPoint1 - data.boxPoint2).ComponentAbs().ComponentDivide(data.resolution);
    }

    TiXmlElement* regularGridElement = dataElement->FirstChildElement("RegularGrid");
    if (regularGridElement)
    {
      // origin
      ParseVector(regularGridElement, "origin" , &data.origin);
      ParseVector(regularGridElement, "spacing", &data.spacing);
    }
  }

  TiXmlElement* meshElement = snapshotElement->FirstChildElement("Mesh");
  if (!meshElement)
  {
    mesh.used = false;
  } else
  {
    ParseString(meshElement, "fileName", &mesh.filename);
  }

  TiXmlElement* contactsElement = snapshotElement->FirstChildElement("Contacts");
  if (!contactsElement)
  {
    contacts.used = false;
  } else
  {
    ParseString(contactsElement, "fileName", &contacts.filename);

    std::string s;
    if (ParseString(contactsElement, "faceTypesToDraw", &s) == TIXML_SUCCESS)
    {
      contacts.faceTypesToDraw = StringToSetOfInt<IndexType>(s);
    }
  }

  TiXmlElement* boundariesElement = snapshotElement->FirstChildElement("Boundaries");
  if (!boundariesElement)
  {
    boundaries.used = false;
  } else
  {
    ParseString(boundariesElement, "fileName", &boundaries.fileName);

    std::string s;
    if (ParseString(boundariesElement, "faceTypesToDraw", &s) == TIXML_SUCCESS)
    {
      boundaries.faceTypesToDraw = StringToSetOfInt<IndexType>(s);
    }
  }

  TiXmlElement* cellInfoElement = snapshotElement->FirstChildElement("CellInfos");
  if (!cellInfoElement)
  {
    cellInfos.used = false;
  } else
  {
    ParseString(cellInfoElement, "fileName", &cellInfos.fileName);
  }
}

template<typename Space>
void DetectorsSettings<Space>::Parse(TiXmlElement* settingsElement)
{
  TiXmlElement* detectorsElement = settingsElement->FirstChildElement("Detectors");
  if (!detectorsElement)
  {
    used = false;
  } else
  {
    ParseString(detectorsElement, "fileName", &filename);
    ParseString(detectorsElement, "locationsFileName", &locationsFileName);
    ParseUnsigned(detectorsElement, "samplesPeriod", &samplesPeriod);
  }
}
