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
#include <cmath>

typedef double MY_TYPE;
#define FORMAT RTAUDIO_FLOAT64
#define PI 3.14159265358979

// ---- Functions from somefunc_2025.cpp ----
int write_buff_dump(MY_TYPE* buff, const int n_buff, MY_TYPE* buff_dump, const int n_buff_dump, int* ind_dump) {
  int i = 0;
  for (i = 0; i < n_buff; i++) {
    if (*ind_dump < n_buff_dump) {
      buff_dump[*ind_dump] = buff[i];
      (*ind_dump)++;
    } else {
      break;
    }
  }
  return i;
}

// Question 12: Autocorrelation function
void autocorrelation(MY_TYPE* x, int Nx, MY_TYPE* y, int Ny) {
  for (int n = 0; n < Ny; n++) {
    y[n] = 0;
    for (int k = n; k < Nx; k++) {
      y[n] += x[k] * x[k - n];
    }
  }
}

// Question 12: Find second maximum for f0 estimation
int find_second_max(MY_TYPE* autocorr, int N, int min_lag, int max_lag) {
  int second_max_index = min_lag;
  MY_TYPE max_value = autocorr[min_lag];

  for (int i = min_lag; i < max_lag && i < N; i++) {
    if (autocorr[i] > max_value) {
      max_value = autocorr[i];
      second_max_index = i;
    }
  }
  return second_max_index;
}

typedef struct mystruct {
  MY_TYPE gain;
  MY_TYPE add;

  // Question 4-6: Audio file loading and playback
  MY_TYPE *fileBuffer;
  long fileLength;
  long fileIndex;

  // Question 7-11: Dump buffers for input and output
  MY_TYPE *dumpBufferInput;
  MY_TYPE *dumpBufferOutput;
  int dumpBufferSize;
  int dumpIndexInput;
  int dumpIndexOutput;

  // Question 12: f0 buffer
  MY_TYPE *dumpBufferF0;
  int dumpBufferF0Size;
  int dumpIndexF0;

  // Sampling rate for f0 calculation
  unsigned int fs;
}* mystruct_t;

void usage( void ) {
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
  std::getline( std::cin, keyHit );
  return i;
}

double streamTimePrintIncrement = 1.0;
double streamTimePrintTime = 1.0;

int inout( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
           double streamTime, RtAudioStreamStatus status, void* data)
{
  if ( status ) std::cout << "Stream over/underflow detected." << std::endl;

  if ( streamTime >= streamTimePrintTime ) {
    std::cout << "streamTime = " << streamTime << std::endl;
    streamTimePrintTime += streamTimePrintIncrement;
  }

  mystruct_t mydata = (mystruct_t) data;
  MY_TYPE* in = (MY_TYPE *)inputBuffer;
  MY_TYPE* out = (MY_TYPE *)outputBuffer;

  // Question 5-6: Load audio file into input buffer (loop playback)
  if (mydata->fileIndex + nBufferFrames < mydata->fileLength) {
    memcpy(in, mydata->fileBuffer + mydata->fileIndex, nBufferFrames * sizeof(MY_TYPE));
    mydata->fileIndex += nBufferFrames;
  }
  else {
    long framesLeft = mydata->fileLength - mydata->fileIndex;
    if (framesLeft > 0) {
      memcpy(in, mydata->fileBuffer + mydata->fileIndex, framesLeft * sizeof(MY_TYPE));
    }
    mydata->fileIndex = 0;
    memcpy(in + framesLeft, mydata->fileBuffer, (nBufferFrames - framesLeft) * sizeof(MY_TYPE));
  }

  // Question 3b: Element-by-element copy (instead of memcpy)
  for (unsigned int i = 0; i < nBufferFrames; i++) {
    out[i] = in[i];
  }

  // Question 8: Write input buffer to dump buffer
  write_buff_dump(in, nBufferFrames, mydata->dumpBufferInput, mydata->dumpBufferSize, &mydata->dumpIndexInput);

  // Question 11: Write output buffer to dump buffer
  write_buff_dump(out, nBufferFrames, mydata->dumpBufferOutput, mydata->dumpBufferSize, &mydata->dumpIndexOutput);

  // Question 12: Autocorrelation and f0 estimation
  MY_TYPE* autocorr = new MY_TYPE[nBufferFrames];
  autocorrelation(in, nBufferFrames, autocorr, nBufferFrames);

  // Find second maximum (f0 estimation)
  // Vocal range: 80-450 Hz -> period range at 44100 Hz: 98-551 samples
  int min_lag = mydata->fs / 450;  // ~98 samples for 44100 Hz
  int max_lag = mydata->fs / 80;   // ~551 samples for 44100 Hz

  int period = find_second_max(autocorr, nBufferFrames, min_lag, max_lag);
  MY_TYPE f0 = (MY_TYPE)mydata->fs / period;

  // Write f0 to dump buffer
  write_buff_dump(&f0, 1, mydata->dumpBufferF0, mydata->dumpBufferF0Size, &mydata->dumpIndexF0);

  delete[] autocorr;

  return 0;
}

int main( int argc, char *argv[] )
{
  unsigned int channels, fs, oDevice = 0, iDevice = 0, iOffset = 0, oOffset = 0;

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

  adac.showWarnings( true );

  unsigned int bufferFrames = 128;
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

  // Question 3c: Using a data structure for userData
  mystruct_t mydata = new mystruct;
  mydata->gain = 0.5;
  mydata->add = 0.1;
  mydata->fs = fs;

  // Question 4: Loading audio file in main()
  const char* audioFile = "../../audio_files/F01_a3_s100_v04.bin";
  FILE *file = fopen(audioFile, "rb");
  if (!file) {
    std::cerr << "Could not open file: " << audioFile << std::endl;
    return 1;
  }
  fseek(file, 0, SEEK_END);
  mydata->fileLength = ftell(file) / sizeof(MY_TYPE);
  fseek(file, 0, SEEK_SET);
  mydata->fileBuffer = new MY_TYPE[mydata->fileLength];
  fread(mydata->fileBuffer, sizeof(MY_TYPE), mydata->fileLength, file);
  mydata->fileIndex = 0;
  fclose(file);

  // Question 7: Initialize dump buffers in main()
  mydata->dumpBufferSize = fs * 6;  // 6 seconds of audio samples
  mydata->dumpBufferInput = new MY_TYPE[mydata->dumpBufferSize];
  mydata->dumpBufferOutput = new MY_TYPE[mydata->dumpBufferSize];

  // f0 buffer: one value per callback, so size = duration * fs / bufferFrames
  mydata->dumpBufferF0Size = 6 * fs / bufferFrames;
  mydata->dumpBufferF0 = new MY_TYPE[mydata->dumpBufferF0Size];

  mydata->dumpIndexInput = 0;
  mydata->dumpIndexOutput = 0;
  mydata->dumpIndexF0 = 0;

  // Question 3d: RTAUDIO_FLOAT64 format (already set)
  if ( adac.openStream( &oParams, &iParams, FORMAT, fs, &bufferFrames, &inout, mydata, &options ) ) {
    goto cleanup;
  }

  if ( adac.isStreamOpen() == false ) goto cleanup;

  std::cout << "\nStream latency = " << adac.getStreamLatency() << " frames" << std::endl;

  if ( adac.startStream() ) goto cleanup;

  char input;
  std::cout << "\nRunning ... press <enter> to quit (buffer frames = " << bufferFrames << ").\n";
  std::cin.get(input);

  if ( adac.isStreamRunning() )
    adac.stopStream();

  // Question 9: Write dump buffers to disk in main()
  {
    FILE *fileInput = fopen("dump_input.bin", "wb");
    fwrite(mydata->dumpBufferInput, sizeof(MY_TYPE), mydata->dumpIndexInput, fileInput);
    fclose(fileInput);

    FILE *fileOutput = fopen("dump_output.bin", "wb");
    fwrite(mydata->dumpBufferOutput, sizeof(MY_TYPE), mydata->dumpIndexOutput, fileOutput);
    fclose(fileOutput);

    // Question 12: Write f0 to disk
    FILE *fileF0 = fopen("dump_f0.bin", "wb");
    fwrite(mydata->dumpBufferF0, sizeof(MY_TYPE), mydata->dumpIndexF0, fileF0);
    fclose(fileF0);

    std::cout << "Wrote " << mydata->dumpIndexInput << " input samples" << std::endl;
    std::cout << "Wrote " << mydata->dumpIndexOutput << " output samples" << std::endl;
    std::cout << "Wrote " << mydata->dumpIndexF0 << " f0 values" << std::endl;
  }

 cleanup:
  if ( adac.isStreamOpen() ) adac.closeStream();

  // Clean up allocated memory
  delete[] mydata->fileBuffer;
  delete[] mydata->dumpBufferInput;
  delete[] mydata->dumpBufferOutput;
  delete[] mydata->dumpBufferF0;
  delete mydata;

  return 0;
}
