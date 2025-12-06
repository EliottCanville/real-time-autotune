/******************************************/
/*
  duplex.cpp
  by Gary P. Scavone, 2006-2019.

  This program opens a duplex stream and passes
  input directly through to the output.
*/
/******************************************/

#include "RtAudio.h"
#include <iostream>
#include <cstdlib>
#include <cstring>

/*
typedef char MY_TYPE;
#define FORMAT RTAUDIO_SINT8

typedef signed short MY_TYPE;
#define FORMAT RTAUDIO_SINT16

typedef S24 MY_TYPE;
#define FORMAT RTAUDIO_SINT24

typedef signed long MY_TYPE;
#define FORMAT RTAUDIO_SINT32

typedef float MY_TYPE;
#define FORMAT RTAUDIO_FLOAT32
*/

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64

typedef struct mystruct {
  MY_TYPE gain;
  MY_TYPE add;

  MY_TYPE *fileBuffer;
  long fileLength;
  long fileIndex;
}* mystruct_t;

void usage( void ) {
  // Error function in case of incorrect command-line
  // argument specifications
  std::cout << "\nuseage: duplex N fs <iDevice> <oDevice> <iChannelOffset> <oChannelOffset>\n";
  std::cout << "    where N = number of channels,\n";
  std::cout << "    fs = the sample rate,\n";
  std::cout << "    iDevice = optional input device index to use (default = 0),\n";
  std::cout << "    oDevice = optional output device index to use (default = 0),\n";
  std::cout << "    iChannelOffset = an optional input channel offset (default = 0),\n";
  std::cout << "    and oChannelOffset = optional output channel offset (default = 0).\n\n";
  exit( 0 );
}

unsigned int getDeviceIndex( std::vector<std::string> deviceNames, bool isInput = false )
{
  unsigned int i;
  std::string keyHit;
  std::cout << '\n';
  for ( i=0; i<deviceNames.size(); i++ )
    std::cout << "  Device #" << i << ": " << deviceNames[i] << '\n';
  do {
    if ( isInput )
      std::cout << "\nChoose an input device #: ";
    else
      std::cout << "\nChoose an output device #: ";
    std::cin >> i;
  } while ( i >= deviceNames.size() );
  std::getline( std::cin, keyHit );  // used to clear out stdin
  return i;
}

double streamTimePrintIncrement = 1.0; // seconds
double streamTimePrintTime = 1.0; // seconds

int inout( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double streamTime, RtAudioStreamStatus status, void* data)
{
  // Since the number of input and output channels is equal, we can do
  // a simple buffer copy operation here.
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  if ( streamTime >= streamTimePrintTime ) {
    std::cout << "streamTime = " << streamTime << std::endl;
    streamTimePrintTime += streamTimePrintIncrement;
  }
  
  mystruct_t mydata = (mystruct_t) data;
  MY_TYPE* in = (MY_TYPE *)inputBuffer;
  MY_TYPE* out = (MY_TYPE *)outputBuffer;

  if (mydata->fileIndex + nBufferFrames < mydata->fileLength) {
    memcpy(in, mydata->fileBuffer + mydata->fileIndex, nBufferFrames * sizeof(MY_TYPE));
    mydata->fileIndex += nBufferFrames;
  }
  else {
    long framesLeft = mydata->fileLength - mydata->fileIndex;
    if (framesLeft > 0) {
      memcpy(in, mydata->fileBuffer + mydata->fileIndex, framesLeft * sizeof(MY_TYPE));
      mydata->fileIndex = 0;
    }
    memcpy(in + framesLeft, mydata->fileBuffer, (nBufferFrames - framesLeft) * sizeof(MY_TYPE));
  }

  for (unsigned int i = 0; i < nBufferFrames; i++)
    // out[i] = mydata->gain * in[i] + mydata->add;
    out[i] = in[i];
  return 0;
}

int main( int argc, char *argv[] )
{
  unsigned int channels, fs, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

  // Minimal command-line checking
  if (argc < 3 || argc > 7 ) usage();

  RtAudio adac;
  std::vector<unsigned int> deviceIds = adac.getDeviceIds();
  if ( deviceIds.size() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 1 );
  }

  channels = (unsigned int) atoi(argv[1]);
  fs = (unsigned int) atoi(argv[2]);
  if ( argc > 3 )
    iDevice = (unsigned int) atoi(argv[3]);
  if ( argc > 4 )
    oDevice = (unsigned int) atoi(argv[4]);
  if ( argc > 5 )
    iOffset = (unsigned int) atoi(argv[5]);
  if ( argc > 6 )
    oOffset = (unsigned int) atoi(argv[6]);

  // Let RtAudio print messages to stderr.
  adac.showWarnings( true );

  // Set the same number of channels for both input and output.
  unsigned int bufferFrames = 512;
  RtAudio::StreamParameters iParams, oParams;
  iParams.nChannels = channels;
  iParams.firstChannel = iOffset;
  oParams.nChannels = channels;
  oParams.firstChannel = oOffset;

  if ( iDevice == 0 )
    iParams.deviceId = adac.getDefaultInputDevice();
  else {
    if ( iDevice >= deviceIds.size() )
      iDevice = getDeviceIndex( adac.getDeviceNames(), true );
    iParams.deviceId = deviceIds[iDevice];
  }
  if ( oDevice == 0 )
    oParams.deviceId = adac.getDefaultOutputDevice();
  else {
    if ( oDevice >= deviceIds.size() )
      oDevice = getDeviceIndex( adac.getDeviceNames() );
    oParams.deviceId = deviceIds[oDevice];
  }
  
  RtAudio::StreamOptions options;
  //options.flags |= RTAUDIO_NONINTERLEAVED;

  mystruct_t mydata = new mystruct;
  mydata->gain = 0.5;
  mydata->add = 0.1;

  // Loading audio file
  FILE *file = fopen("/Users/eliott/Documents/ETUDES/ENSE3-PHELMA/3A/S9/Temps Réel/BE_tstr_2025_v1/audio_files/F01_a3_s100_v04.bin", "rb");
  fseek(file, 0, SEEK_END);
  mydata->fileLength = ftell(file) / sizeof(MY_TYPE);
  fseek(file, 0, SEEK_SET);
  mydata->fileBuffer = new MY_TYPE[mydata->fileLength];
  size_t result = fread(mydata->fileBuffer, sizeof(MY_TYPE), mydata->fileLength, file);
  mydata->fileIndex = 0;
  fclose(file);


  if ( adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, mydata, &options ) ) {
    goto cleanup;
  }

  if ( adac.isStreamOpen() == false ) goto cleanup;

  // Test RtAudio functionality for reporting latency.
  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  if ( adac.startStream() ) goto cleanup;

  char input;
  std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
  std::cin.get(input);

  // Stop the stream.
  if ( adac.isStreamRunning() )
    adac.stopStream();

 cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  return 0;
}
