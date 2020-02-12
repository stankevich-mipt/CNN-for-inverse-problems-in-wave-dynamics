#pragma once
#include "ParserUtil.h"

template<typename Space>
struct TaskSettings
{
  SPACE_TYPEDEFS

  TaskSettings():
    destinationTime(-1.0)
  {
  }

  struct BoxStateInfo
  {
    BoxStateInfo():
      boxPoint1(Vector::zeroVector()),
      boxPoint2(Vector::zeroVector()),
      velocity(Vector::zeroVector())
    {
    }
    Vector boxPoint1;
    Vector boxPoint2;
    Vector velocity;
  };
  std::vector<BoxStateInfo> boxStateInfos;

  struct BerlageStateInfo
  {
    BerlageStateInfo():
      point(Vector::zeroVector()),
      waveVector(Vector::xAxis()),
      waveLength(1.0),
      shear(false)
    {
    }
    Vector point;
    Vector waveVector;
    Scalar waveLength;
    bool shear;
  };
  std::vector<BerlageStateInfo> berlageStateInfos;

  struct RickerStateInfo
  {
    RickerStateInfo():
      point(Vector::zeroVector()),
      waveVector(Vector::xAxis()),
      waveLength(1.0)
    {
    }
    Vector point;
    Vector waveVector;
    Scalar waveLength;
  };
  std::vector<RickerStateInfo> rickerStateInfos;

  struct AcousticPlaneWaveStateInfo
  {
    AcousticPlaneWaveStateInfo():
      point(Vector::zeroVector()),
      waveVector(Vector::xAxis()),
      waveLength(1.0)
    {
    }
    Vector point;
    Vector waveVector;
    Scalar waveLength;
  };
  std::vector<AcousticPlaneWaveStateInfo> acousticPlaneWaveStateInfos;

  struct BallStateInfo
  {
    BallStateInfo():
      point(Vector::zeroVector()),
      radialVelocity(1.0),
      pressure(1.0),
      waveLength(1.0)
    {
    }
    Vector point;
    Scalar radialVelocity;
    Scalar pressure;
    Scalar waveLength;
  };
  std::vector<BallStateInfo> ballStateInfos;

  struct RadialStateInfo
  {
    RadialStateInfo():
      point(Vector::zeroVector()),
      radialVelocity(1.0),
      waveLength(1.0)
    {
    }
    Vector point;
    Scalar radialVelocity;
    Scalar waveLength;
  };
  std::vector<RadialStateInfo> radialStateInfos;

  struct BagelStateInfo
  {
    BagelStateInfo():
      r1(1.0),
      r2(2.0),
      point(Vector::zeroVector()),
      magnitude(1.0)
    {
    }
    Scalar r1;
    Scalar r2;
    Vector point;
    Scalar magnitude;
  };
  std::vector<BagelStateInfo> bagelStateInfos;

  struct DirectionalExplosionStateInfo
  {
    DirectionalExplosionStateInfo():
      point(Vector::zeroVector()),
      velocity(Vector::xAxis()),
      waveLength(1.0)
    {
    }
    Vector point;
    Vector velocity;
    Scalar waveLength;
  };
  std::vector<DirectionalExplosionStateInfo> directionalExplosionStateInfos;

 struct MonochromaticPlaneWaveStateInfo
  {
    MonochromaticPlaneWaveStateInfo():
      normal(Vector::xAxis()),
      waveLength(1.0),
      initialPhase(0.0),
      shear(false)
    {
    }
    Vector normal;
    Scalar waveLength;
    Scalar initialPhase;
    bool shear;
  };
  std::vector<MonochromaticPlaneWaveStateInfo> monochromaticPlaneWaveStateInfos;

  struct ConstantPlaneWaveStateInfo
  {
    ConstantPlaneWaveStateInfo():
      velocity(Vector::xAxis()),
      center(Vector::zeroVector()),
      waveLength(1.0),
      shear(false)
    {
    }
    Vector velocity;
    Vector center;
    Scalar waveLength;
    bool shear;
  };
  std::vector<ConstantPlaneWaveStateInfo> constantPlaneWaveStateInfos;

  struct IniState
  {
    enum IniStateTypes
    {
      Box,
      Berlage,
      Ricker,
      Radial,
      Ball,
      Bagel,
      DirectionalExplosion,
      MonochromaticPlaneWave,
      ConstantPlaneWave,
      AcousticPlaneWave,
      Undefined
    };

    IniStateTypes type;
    IndexType infoIndex;
  };
  std::vector<IniState> iniStates;

  // source terms
  struct ConstantAccelerationSourceInfo
  {
    ConstantAccelerationSourceInfo():
      acceleration(Vector::zeroVector())
    {}
    Vector acceleration;
  };
  std::vector<ConstantAccelerationSourceInfo> constantAccelerationSourceInfos;

  struct SourceTerm
  {
    enum SourceTermsTypes
    {
      ConstantAcceleration,
      Undefined
    };

    SourceTermsTypes type;
    IndexType infoIndex;
  };
  std::vector<SourceTerm> sourceTerms;

  // point sources
  struct ForcePointSourceInfo
  {
    ForcePointSourceInfo():
      point(Vector::zeroVector()),
      peakFrequency(Scalar(1.0)), 
      acceleration(Vector::xAxis()),
      latency(Scalar(0.0))
    {}
    Vector point;
    Scalar peakFrequency;
    Vector acceleration;
    Scalar latency;
  };
  std::vector<ForcePointSourceInfo> forcePointSourceInfos;

  struct MonopoleSourceInfo
  {
    MonopoleSourceInfo():
      point(Vector::zeroVector()),
      peakFrequency(Scalar(1.0)),
      pressure(Scalar(1.0)),
      latency(Scalar(0.0))
    {}
    Vector point;
    Scalar peakFrequency;
    Scalar pressure;
    Scalar latency;
  };
  std::vector<MonopoleSourceInfo> monopoleSourceInfos;

  struct PointSource
  {
    enum PointSourcesTypes
    {
      Force,
      Monopole,
      Undefined
    };
    PointSourcesTypes type;
    IndexType infoIndex;
  };
  std::vector<PointSource> pointSources;

  Scalar destinationTime;

  void Parse(TiXmlElement *taskElement);
};

template<typename Space>
void TaskSettings<Space>::Parse(TiXmlElement *taskElement)
{
  ParseScalar(taskElement, "destinationTime", &destinationTime);

  TiXmlElement *iniStateElement = taskElement->FirstChildElement("IniState");

  if(iniStateElement)
  {
    TiXmlElement *boxStateElement = iniStateElement->FirstChildElement("BoxState");
    while(boxStateElement)
    {
      IniState state;
      state.type = IniState::Box;
      state.infoIndex = boxStateInfos.size();

      BoxStateInfo boxInfo;
      ParseVector(boxStateElement, "boxPoint1",  &(boxInfo.boxPoint1));
      ParseVector(boxStateElement, "boxPoint2",  &(boxInfo.boxPoint2));
      ParseVector(boxStateElement, "velocity",   &(boxInfo.velocity ));
      boxStateInfos.push_back(boxInfo);

      iniStates.push_back(state);

      boxStateElement = boxStateElement->NextSiblingElement("BoxState");
    }

    TiXmlElement *rickerStateElement = iniStateElement->FirstChildElement("RickerWave");
    while(rickerStateElement)
    {
      IniState state;
      state.type = IniState::Ricker;
      state.infoIndex = rickerStateInfos.size();

      RickerStateInfo rickerInfo;
      ParseVector(rickerStateElement, "point",       &(rickerInfo.point      ));
      ParseVector(rickerStateElement, "waveVector",  &(rickerInfo.waveVector ));
      ParseScalar(rickerStateElement, "waveLength",  &(rickerInfo.waveLength ));
      rickerStateInfos.push_back(rickerInfo);

      iniStates.push_back(state);

      rickerStateElement = rickerStateElement->NextSiblingElement("RickerWave");
    }

    TiXmlElement *berlageStateElement = iniStateElement->FirstChildElement("BerlageWave");
    while(berlageStateElement)
    {
      IniState state;
      state.type = IniState::Berlage;
      state.infoIndex = berlageStateInfos.size();

      BerlageStateInfo berlageInfo;
      ParseVector(berlageStateElement, "point",       &(berlageInfo.point      ));
      ParseVector(berlageStateElement, "waveVector",  &(berlageInfo.waveVector ));
      ParseScalar(berlageStateElement, "waveLength",  &(berlageInfo.waveLength ));
      ParseBool  (berlageStateElement, "shear",       &(berlageInfo.shear      ));
      berlageStateInfos.push_back(berlageInfo);

      iniStates.push_back(state);

      berlageStateElement = berlageStateElement->NextSiblingElement("BerlageWave");
    }

    TiXmlElement *acousticPlaneWaveStateElement = iniStateElement->FirstChildElement("AcousticPlaneWave");
    while (acousticPlaneWaveStateElement)
    {
      IniState state;
      state.type = IniState::AcousticPlaneWave;
      state.infoIndex = acousticPlaneWaveStateInfos.size();

      AcousticPlaneWaveStateInfo acousticPlaneWaveInfo;
      ParseVector(acousticPlaneWaveStateElement, "point",      &(acousticPlaneWaveInfo.point));
      ParseVector(acousticPlaneWaveStateElement, "waveVector", &(acousticPlaneWaveInfo.waveVector));
      ParseScalar(acousticPlaneWaveStateElement, "waveLength", &(acousticPlaneWaveInfo.waveLength));

      acousticPlaneWaveStateInfos.push_back(acousticPlaneWaveInfo);
      iniStates.push_back(state);
      acousticPlaneWaveStateElement = acousticPlaneWaveStateElement->NextSiblingElement("AcousticPlaneWave");
    }

    TiXmlElement *constantPlaneWaveStateElement = iniStateElement->FirstChildElement("ConstantPlaneWave");
    while (constantPlaneWaveStateElement)
    {
      IniState state;
      state.type = IniState::ConstantPlaneWave;
      state.infoIndex = constantPlaneWaveStateInfos.size();

      ConstantPlaneWaveStateInfo constantPlaneWaveInfo;
      ParseVector(constantPlaneWaveStateElement, "velocity", &(constantPlaneWaveInfo.velocity));
      ParseVector(constantPlaneWaveStateElement, "center", &(constantPlaneWaveInfo.center));
      ParseScalar(constantPlaneWaveStateElement, "waveLength", &(constantPlaneWaveInfo.waveLength));
      ParseBool(constantPlaneWaveStateElement, "shear", &(constantPlaneWaveInfo.shear));

      constantPlaneWaveStateInfos.push_back(constantPlaneWaveInfo);
      iniStates.push_back(state);
      constantPlaneWaveStateElement = constantPlaneWaveStateElement->NextSiblingElement("ConstantPlaneWave");
    }

    TiXmlElement *radialWaveElement = iniStateElement->FirstChildElement("RadialWave");
    while(radialWaveElement)
    {
      IniState state;
      state.type = IniState::Radial;
      state.infoIndex = radialStateInfos.size();

      RadialStateInfo radialInfo;
      ParseVector(radialWaveElement, "point",  &(radialInfo.point));
      ParseScalar(radialWaveElement,  "radialVelocity",  &(radialInfo.radialVelocity));
      ParseScalar(radialWaveElement,  "waveLength",   &(radialInfo.waveLength));
      radialStateInfos.push_back(radialInfo);

      iniStates.push_back(state);

      radialWaveElement = radialWaveElement->NextSiblingElement("RadialWave");
    }

    TiXmlElement *ballWaveElement = iniStateElement->FirstChildElement("BallWave");
    while(ballWaveElement)
    {
      IniState state;
      state.type = IniState::Ball;
      state.infoIndex = ballStateInfos.size();

      BallStateInfo ballInfo;
      ParseVector(ballWaveElement, "point",  &(ballInfo.point));
      ParseScalar(ballWaveElement,  "radialVelocity",  &(ballInfo.radialVelocity));
      ParseScalar(ballWaveElement,  "waveLength",   &(ballInfo.waveLength));
      ParseScalar(ballWaveElement,  "pressure",   &(ballInfo.pressure));
      ballStateInfos.push_back(ballInfo);

      iniStates.push_back(state);

      ballWaveElement = ballWaveElement->NextSiblingElement("BallWave");
    }

    TiXmlElement *bagelWaveElement = iniStateElement->FirstChildElement("BagelWave");
    while(bagelWaveElement)
    {
      IniState state;
      state.type = IniState::Bagel;
      state.infoIndex = bagelStateInfos.size();

      BagelStateInfo bagelInfo;
      ParseVector(bagelWaveElement, "point",  &(bagelInfo.point));
      ParseScalar(bagelWaveElement,  "magnitude",  &(bagelInfo.magnitude));
      ParseScalar(bagelWaveElement,  "r1",   &(bagelInfo.r1));
      ParseScalar(bagelWaveElement,  "r2",   &(bagelInfo.r2));
      bagelStateInfos.push_back(bagelInfo);

      iniStates.push_back(state);

      bagelWaveElement = bagelWaveElement->NextSiblingElement("BagelWave");
    }

    TiXmlElement *directionalExplosionElement = iniStateElement->FirstChildElement("DirectionalExplosion");
    while(directionalExplosionElement)
    {
      IniState state;
      state.type = IniState::DirectionalExplosion;
      state.infoIndex = directionalExplosionStateInfos.size();

      DirectionalExplosionStateInfo explosionInfo;
      ParseVector(directionalExplosionElement, "point",      &(explosionInfo.point));
      ParseVector(directionalExplosionElement, "velocity",   &(explosionInfo.velocity));
      ParseScalar(directionalExplosionElement, "waveLength", &(explosionInfo.waveLength));
      directionalExplosionStateInfos.push_back(explosionInfo);

      iniStates.push_back(state);

      directionalExplosionElement = directionalExplosionElement->NextSiblingElement("DirectionalExplosion");
    }
  }

  TiXmlElement *sourceTermsElement = taskElement->FirstChildElement("SourceTerms");

  if(sourceTermsElement)
  {
    TiXmlElement *constantAccelerationElement = sourceTermsElement->FirstChildElement("ConstantAcceleration");
    while(constantAccelerationElement)
    {
      SourceTerm terms;
      terms.type = SourceTerm::ConstantAcceleration;
      terms.infoIndex = constantAccelerationSourceInfos.size();

      ConstantAccelerationSourceInfo info;
      ParseVector(constantAccelerationElement, "acceleration", &(info.acceleration));
      constantAccelerationSourceInfos.push_back(info);

      sourceTerms.push_back(terms);

      constantAccelerationElement = constantAccelerationElement->NextSiblingElement("ConstantAcceleration");
    }
  }

  TiXmlElement* pointSourcesElement = taskElement->FirstChildElement("PointSources");

  if (pointSourcesElement)
  {
    TiXmlElement* forcePointSourceElement = pointSourcesElement->FirstChildElement("ForceSource");
    while (forcePointSourceElement)
    {
      PointSource source;
      source.type = PointSource::Force;
      source.infoIndex = forcePointSourceInfos.size();

      ForcePointSourceInfo info;
      ParseVector(forcePointSourceElement, "point",         &(info.point));
      ParseVector(forcePointSourceElement, "acceleration",  &(info.acceleration));
      ParseScalar(forcePointSourceElement, "peakFrequency", &(info.peakFrequency));
      ParseScalar(forcePointSourceElement, "latency",       &(info.latency));
      forcePointSourceInfos.push_back(info);

      pointSources.push_back(source);

      forcePointSourceElement = forcePointSourceElement->NextSiblingElement("ForceSource");
    }

    TiXmlElement* monopoleSourceElement = pointSourcesElement->FirstChildElement("MonopoleSource");
    while (monopoleSourceElement)
    {
      PointSource source;
      source.type = PointSource::Monopole;
      source.infoIndex = monopoleSourceInfos.size();

      MonopoleSourceInfo info;
      ParseVector(monopoleSourceElement, "point", &(info.point));
      ParseScalar(monopoleSourceElement, "pressure", &(info.pressure));
      ParseScalar(monopoleSourceElement, "peakFrequency", &(info.peakFrequency));
      ParseScalar(monopoleSourceElement, "latency", &(info.latency));
      monopoleSourceInfos.push_back(info);

      pointSources.push_back(source);

      monopoleSourceElement = monopoleSourceElement->NextSiblingElement("MonopoleSource");
    }
  }
}
