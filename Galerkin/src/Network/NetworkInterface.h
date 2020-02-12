#include "mpi.h"
#include <vector>
struct NotifyListener
{
  virtual ~NotifyListener() = default;
  virtual void OnNotify() = 0;
};

struct ReceiveListener
{
  virtual ~ReceiveListener() = default;
  virtual void OnDataReceive(void *data, int size, int hostId) = 0;
};

class NetworkInterface
{
public:
  NetworkInterface();
  void init(int *i_argc, char*** i_argv);
  void finalize();
  int getID();
  double getWTime();
  int getHostsCount();
  bool checkAllReadyAsync();
  void barrier();
  void setReceiveListener(ReceiveListener* io_receiveListener);
  void notifyReady();
  void sendDataSync(void* i_buffer, int i_bufSize, int i_hostId, NotifyListener* i_notifyListener);
  void sendDataAsync(void* i_buffer, int i_bufSize, int i_hostId, NotifyListener* i_notifyListener);
  void update(/*int* receivedNotifications*/);
  
  void waitAll();
  void initializeSync();
  
  ReceiveListener* const GetReceiveListener() const
  {
    return receiveListener;
  }

  struct DataNotifyListener: public NotifyListener
  {
    DataNotifyListener(): 
      NotifyListener(),
      notificationsRemained(0)
    {}
  
    virtual void OnNotify() override
    {
      notificationsRemained--;
    }
    int notificationsRemained;
  };

  template <typename DataType, typename Comparator>
  struct DataReceiveListener: public ReceiveListener
  {
    DataReceiveListener(DataType data, Comparator comparator): 
      ReceiveListener(),
      commonData(data),
      comparator(comparator)
    {}
    virtual void OnDataReceive(void *data, int size, int hostId) override
    {
      DataType value = *((DataType*)data);
      commonData = comparator(value, commonData);
    }
    DataType commonData;
    Comparator comparator;
  };

  template <typename DataType, typename Comparator>
  DataType Negotiate(DataType data, const Comparator comparator)
  {
    initializeSync();
    DataNotifyListener stepNotifyListener;
    DataReceiveListener<DataType, Comparator> stepReceiveListener(data, comparator);

    stepReceiveListener.commonData = data;
    ReceiveListener* const previousReceiveListener = this->GetReceiveListener();
    this->setReceiveListener(&stepReceiveListener);

    DataType *sendBuf = new DataType[this->getHostsCount()];
    for(int dstHost = 0; dstHost < this->getHostsCount(); dstHost++)
    {
      sendBuf[dstHost] = data;
      if(dstHost != this->getID())
      {
        this->sendDataAsync(sendBuf + dstHost, sizeof(DataType), int(dstHost), &stepNotifyListener);
        stepNotifyListener.notificationsRemained++;
      }
    }

    // printf("Node %d is waiting for nodes to receive its message\n", this->getID());
    bool notified = false;
    while(!this->checkAllReadyAsync())
    {
      this->update();
      if(stepNotifyListener.notificationsRemained == 0 && !notified)
      {
        // printf("Node %d has received all notifications\n", this->getID());
        this->notifyReady();
        notified = true;
      }
    }
    // printf("Node %d is at barrier\n", this->getID());
    this->barrier();

    delete [] sendBuf;
    this->setReceiveListener(previousReceiveListener);
    waitAll();
    return stepReceiveListener.commonData;
  }

private:
  enum MessageTags
  {
    ReadyNotification,
    ReceivedNotification,
    Data
  };
  void sendReceivedNotification(int i_dest);
  bool blokingReceive(int &infoType, int &sourceId);
  bool nonBlokingReceive(int &infoType, int &sourceId);
private:
  //process info
  int rank;
  //for handling receive event
  ReceiveListener* receiveListener;
  NotifyListener* notifyListener;
  //buffers
  std::vector<char> receiveBuffer;
  int notifyBuffer;//type indicator of notification
  int currentBufSize;//size of info in receiveBuffer
  int infoType;
  //information for handling interaction
  int lastSource;//id of last process netInerface interact with
  MPI_Status status;
  int numReadyNotifications;//number of received ready notifications
  int numAllProc;
  bool readyForBarrier;

  std::vector<MPI_Request> requests;
  int requestsCount;

  MPI_Request* GetRequest()
  {
    requestsCount++;
    if (size_t(requestsCount) > requests.size()) requests.resize(requestsCount);
    requests[requestsCount - 1] = MPI_REQUEST_NULL;
    return &requests[requestsCount - 1];
  }
};
