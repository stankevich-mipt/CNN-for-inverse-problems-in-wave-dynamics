#include <iostream>
#include <assert.h>

#include "NetworkInterface.h"
NetworkInterface::NetworkInterface()
{
  numReadyNotifications = 0;
  notifyBuffer = 0;
  rank = -1;
  readyForBarrier = false;
  receiveListener = nullptr;
  notifyListener  = nullptr;
}
/*
denote the start of parallel program
*/
void NetworkInterface::init(int *i_argc, char*** i_argv)
{
  MPI_Init(i_argc, i_argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numAllProc);
}
/*
denote the end of parallel program
*/
void NetworkInterface::finalize()
{
  MPI_Finalize();
}
/*
return identification number in according to MPI mapping
*/
int NetworkInterface::getID()
{
  if(rank < 0)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  return rank;
}

int NetworkInterface::getHostsCount()
{
  return numAllProc;
}

double NetworkInterface::getWTime()
{
  return MPI_Wtime();
}
/*
check whether current host received all READY_NOTIFICATION messages
*/
bool NetworkInterface::checkAllReadyAsync()
{
  if (numReadyNotifications == numAllProc)
  {
    numReadyNotifications = 0;
    return true;
  }
  return false;
}
void NetworkInterface::barrier()
{
  MPI_Barrier(MPI_COMM_WORLD);
}
void NetworkInterface::setReceiveListener(ReceiveListener* io_receiveListener)
{
  this->receiveListener = io_receiveListener;
}
/*
we will send notification(READY_NOTIFICATION) that have already received all notifications that we need
to all workers in the network
( it isn't efficient way, some modifications can be done to increase efficiency )
*/
void NetworkInterface::notifyReady()
{
  for (int j = 0; j < numAllProc; j++)
  {
    MPI_Isend(&notifyBuffer, 1, MPI_INT, j, ReadyNotification, MPI_COMM_WORLD, GetRequest());
  }
}
/*
return: always true
*/
bool NetworkInterface::blokingReceive(int &infoType, int &sourceId)
{
  int count;
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);//tag defines the process in receiving and sending( rank = tag )
  infoType = status.MPI_TAG;
  switch(status.MPI_TAG)
  {
    case ReceivedNotification:
    {
      MPI_Get_count(&status, MPI_INT, &count);
      MPI_Recv(&notifyBuffer, count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      currentBufSize = 0;
    }break;

    case ReadyNotification:
    {
      MPI_Get_count(&status, MPI_INT, &count);
      MPI_Recv(&notifyBuffer, count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      currentBufSize = 0;
    }break;

    case Data:
    {
      MPI_Get_count(&status, MPI_CHAR, &count);
      if(receiveBuffer.size() < size_t(count)) receiveBuffer.resize(count);
      MPI_Recv(&(receiveBuffer[0]), count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
      currentBufSize = count;
    }break;
  }
  sourceId = status.MPI_SOURCE;
  return true;
}
/*
return: true  if we actually receive some info
    false if we there isn't information to receive
*/
bool NetworkInterface::nonBlokingReceive(int &infoType, int &sourceId)
{
  int flag;
  int count;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);//tag defines the process in receiving and sending

  if(flag)
  {
    infoType = status.MPI_TAG;
    sourceId = status.MPI_SOURCE;

    switch(status.MPI_TAG)
    {
      case ReceivedNotification:
      {
        MPI_Get_count(&status, MPI_INT, &count);
        assert(count == 1);
        MPI_Recv(&notifyBuffer, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }break;

      case ReadyNotification:
      {
        MPI_Get_count(&status, MPI_INT, &count);
        assert(count == 1);
        MPI_Recv(&notifyBuffer, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }break;

      case Data:
      {
        MPI_Get_count(&status, MPI_CHAR, &count);
        if(receiveBuffer.size() < size_t(count)) receiveBuffer.resize(count * sizeof(char));
        MPI_Recv(&(receiveBuffer[0]), count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        currentBufSize = count;
      }break;
    }
  }
  return flag ? true : false;
}

void NetworkInterface::sendDataSync(void *i_buffer, int i_bufSize, int i_hostId, NotifyListener* i_notifyListener)
{
  this->notifyListener = i_notifyListener;
  MPI_Send(i_buffer, i_bufSize, MPI_CHAR, i_hostId, Data, MPI_COMM_WORLD);
}

void NetworkInterface::sendDataAsync(void *i_buffer, int i_bufSize, int i_hostId, NotifyListener* i_notifyListener)
{
  this->notifyListener = i_notifyListener;
  MPI_Isend(i_buffer, i_bufSize, MPI_CHAR, i_hostId, Data, MPI_COMM_WORLD, GetRequest());
}

void NetworkInterface::sendReceivedNotification(int i_dest)
{
  MPI_Isend(&notifyBuffer, 1, MPI_INT, i_dest, ReceivedNotification, MPI_COMM_WORLD, GetRequest());
}
/*
it is the main method for interacting with other hosts
this method receives three types of messages and handle it in appropriate way
1)MESSAGE - message for updating information in region
2)SIMPLE_NOTIFICATION - message for denoting that we receive MESSAGE
3)READY_NOTIFICATION - message for denoting that we have already received all SIMPLE_NOTIFICATION
*/
void NetworkInterface::update(/*int* receivedNotifications*/)
{
//  if(blokingReceive())//we get new information
  int infoType;
  int sourceId;
  if(nonBlokingReceive(infoType, sourceId))//we get new information
  {
    switch(infoType)
    {
      case Data:
        sendReceivedNotification(sourceId);
        receiveListener->OnDataReceive(&(receiveBuffer[0]), currentBufSize, sourceId);
        break;
      case ReceivedNotification:
        notifyListener->OnNotify();
        break;
      case ReadyNotification:
        numReadyNotifications++;
        break;
    }
  }
}

void NetworkInterface::waitAll()
{
  if (requestsCount > 0)
  {
    MPI_Waitall(requestsCount, &(requests[0]), MPI_STATUSES_IGNORE);
  }
}

void NetworkInterface::initializeSync()
{
  requestsCount = 0;
}

