template<typename Space>
struct ScheduleSettings
{
  SPACE_TYPEDEFS

  ScheduleSettings(): 
    domainsCount(1)
  {}
  IndexType domainsCount;

  struct NodeSchedule
  {
    std::vector<IndexType> domainsIndices;
  };

  // lists of domains which will be computed on the i-th computation node
  std::vector<NodeSchedule> nodesSchedule;


  void Parse(TiXmlElement *scheduleElement);
};

template<typename Space>
void ScheduleSettings<Space>::Parse(TiXmlElement *scheduleElement)
{
  if (scheduleElement->QueryUnsignedAttribute("domainsCount", reinterpret_cast<unsigned int*>(&domainsCount)) != TIXML_SUCCESS)
  {
    std::cerr << "There is no Schedule.domainsCount attribute";
    throw;
  }
   
  TiXmlElement* nodeElement = scheduleElement->FirstChildElement("NodeSchedule");

  if (!nodeElement)
  {
    // i-th node <-> i-th domain
    IndexType nodesCount = domainsCount;
    nodesSchedule.resize(nodesCount);
    for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
    {
      IndexType nodeIndex = domainIndex;
      nodesSchedule[nodeIndex].domainsIndices.push_back(domainIndex);
    }
  }

  while (nodeElement)
  {
    unsigned int nodeIndex;

    if (nodeElement->QueryUnsignedAttribute("node", static_cast<unsigned int*>(&nodeIndex)) != TIXML_SUCCESS)
    {
      std::cerr << "There is no attribute NodeSchedule.node";
      throw;
    }
    nodesSchedule.resize(std::max(static_cast<size_t>(nodeIndex + 1), nodesSchedule.size()));
    std::string domains;

    if (nodeElement->QueryStringAttribute("domains", &domains) != TIXML_SUCCESS)
    {
      std::cerr << "There is no attribute NodeSchedule.domains";
      throw;
    }

    if (domains == "-1")
    {
      for (IndexType domainIndex = 0; domainIndex < domainsCount; ++domainIndex)
      {
        nodesSchedule[nodeIndex].domainsIndices.push_back(domainIndex);
      }
    } else
    {
      std::stringstream s(domains);
      IndexType domainIndex;
      while (!s.eof())
      {
        s >> domainIndex;
        nodesSchedule[nodeIndex].domainsIndices.push_back(domainIndex);
      }
    }
    nodeElement = nodeElement->NextSiblingElement("NodeSchedule");
  }
}
