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

// ---- FFT-related functions (from somefunc_2025.cpp) ----
MY_TYPE *_sintbl = 0;
int maxfftsize = 0;

char *getmem(int leng, unsigned size) {
  char *p = NULL;
  if ((p = (char *)calloc(leng, size)) == NULL){
    fprintf(stderr, "Memory allocation error !\n");
    exit(3);
  }
  return (p);
}

MY_TYPE *dgetmem(int leng) {
  return ((MY_TYPE *)getmem(leng, sizeof(MY_TYPE)));
}

int get_nextpow2(int n) {
  int k = 1;
  while (k < n){
    k *= 2;
  }
  return k;
}

static int checkm(const int m) {
  int k;
  for (k = 4; k <= m; k <<= 1) {
    if (k == m)
      return (0);
  }
  fprintf(stderr, "fft : m must be a integer of power of 2! (m=%i)\n",m);
  return (-1);
}

int fft(MY_TYPE *x, MY_TYPE *y, const int m) {
  int j, lmx, li;
  MY_TYPE *xp, *yp;
  MY_TYPE *sinp, *cosp;
  int lf, lix, tblsize;
  int mv2, mm1;
  MY_TYPE t1, t2;
  MY_TYPE arg;

  if (checkm(m))
    return (-1);

  if ((_sintbl == 0) || (maxfftsize < m)) {
    tblsize = m - m / 4 + 1;
    arg = PI / m * 2;
    if (_sintbl != 0)
      free(_sintbl);
    _sintbl = sinp = dgetmem(tblsize);
    *sinp++ = 0;
    for (j = 1; j < tblsize; j++)
      *sinp++ = sin(arg * (MY_TYPE) j);
    _sintbl[m / 2] = 0;
    maxfftsize = m;
  }

  lf = maxfftsize / m;
  lmx = m;

  for (;;) {
    lix = lmx;
    lmx /= 2;
    if (lmx <= 1)
      break;
    sinp = _sintbl;
    cosp = _sintbl + maxfftsize / 4;
    for (j = 0; j < lmx; j++) {
      xp = &x[j];
      yp = &y[j];
      for (li = lix; li <= m; li += lix) {
        t1 = *(xp) - *(xp + lmx);
        t2 = *(yp) - *(yp + lmx);
        *(xp) += *(xp + lmx);
        *(yp) += *(yp + lmx);
        *(xp + lmx) = *cosp * t1 + *sinp * t2;
        *(yp + lmx) = *cosp * t2 - *sinp * t1;
        xp += lix;
        yp += lix;
      }
      sinp += lf;
      cosp += lf;
    }
    lf += lf;
  }

  xp = x;
  yp = y;
  for (li = m / 2; li--; xp += 2, yp += 2) {
    t1 = *(xp) - *(xp + 1);
    t2 = *(yp) - *(yp + 1);
    *(xp) += *(xp + 1);
    *(yp) += *(yp + 1);
    *(xp + 1) = t1;
    *(yp + 1) = t2;
  }

  j = 0;
  xp = x;
  yp = y;
  mv2 = m / 2;
  mm1 = m - 1;
  for (lmx = 0; lmx < mm1; lmx++) {
    if ((li = lmx - j) < 0) {
      t1 = *(xp);
      t2 = *(yp);
      *(xp) = *(xp + li);
      *(yp) = *(yp + li);
      *(xp + li) = t1;
      *(yp + li) = t2;
    }
    li = mv2;
    while (li <= j) {
      j -= li;
      li /= 2;
    }
    j += li;
    xp = x + j;
    yp = y + j;
  }

  return (0);
}

int fftr(MY_TYPE *x, MY_TYPE *y, const int m) {
  int i, j;
  MY_TYPE *xp, *yp, *xq;
  MY_TYPE *yq;
  int mv2, n, tblsize;
  MY_TYPE xt, yt, *sinp, *cosp;
  MY_TYPE arg;

  mv2 = m / 2;

  xq = xp = x;
  yp = y;
  for (i = mv2; --i >= 0;) {
    *xp++ = *xq++;
    *yp++ = *xq++;
  }

  if (fft(x, y, mv2) == -1)
    return (-1);

  if ((_sintbl == 0) || (maxfftsize < m)) {
    tblsize = m - m / 4 + 1;
    arg = PI / m * 2;
    if (_sintbl != 0)
      free(_sintbl);
    _sintbl = sinp = dgetmem(tblsize);
    *sinp++ = 0;
    for (j = 1; j < tblsize; j++)
      *sinp++ = sin(arg * (MY_TYPE) j);
    _sintbl[m / 2] = 0;
    maxfftsize = m;
  }

  n = maxfftsize / m;
  sinp = _sintbl;
  cosp = _sintbl + maxfftsize / 4;

  xp = x;
  yp = y;
  xq = xp + m;
  yq = yp + m;
  *(xp + mv2) = *xp - *yp;
  *xp = *xp + *yp;
  *(yp + mv2) = *yp = 0;

  for (i = mv2, j = mv2 - 2; --i; j -= 2) {
    ++xp;
    ++yp;
    sinp += n;
    cosp += n;
    yt = *yp + *(yp + j);
    xt = *xp - *(xp + j);
    *(--xq) = (*xp + *(xp + j) + *cosp * yt - *sinp * xt) * 0.5;
    *(--yq) = (*(yp + j) - *yp + *sinp * yt + *cosp * xt) * 0.5;
  }

  xp = x + 1;
  yp = y + 1;
  xq = x + m;
  yq = y + m;

  for (i = mv2; --i;) {
    *xp++ = *(--xq);
    *yp++ = -(*(--yq));
  }

  return (0);
}

static MY_TYPE *hanning(MY_TYPE *w, const int leng) {
  int i;
  MY_TYPE arg;
  MY_TYPE *p;

  arg = 2*PI / (leng - 1);
  for (p = w, i = 0; i < leng; i++)
    *p++ = 0.5 * (1 - cos(i * arg));

  return (w);
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

  // Question 15-18: FFT and harmonic analysis
  int n_fft;              // FFT size
  int max_harmonics;      // Maximum number of harmonics to track
  MY_TYPE *phases;        // Phase accumulator for each harmonic (for Q18)
  long sampleCounter;     // Sample counter for phase continuity
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

  // Question 8: Write input buffer to dump buffer
  write_buff_dump(in, nBufferFrames, mydata->dumpBufferInput, mydata->dumpBufferSize, &mydata->dumpIndexInput);

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

  // Question 15-16: FFT and harmonic analysis
  int n_fft = mydata->n_fft;

  // Allocate FFT buffers (real and imaginary parts)
  MY_TYPE* fft_real = new MY_TYPE[n_fft];
  MY_TYPE* fft_imag = new MY_TYPE[n_fft];

  // Copy input to FFT buffer and zero-pad if necessary
  for (unsigned int i = 0; i < nBufferFrames && i < (unsigned int)n_fft; i++) {
    fft_real[i] = in[i];
    fft_imag[i] = 0.0;
  }
  for (int i = nBufferFrames; i < n_fft; i++) {
    fft_real[i] = 0.0;
    fft_imag[i] = 0.0;
  }

  // Calculate FFT
  fftr(fft_real, fft_imag, n_fft);

  // Question 15: Calculate DFT bin closest to f0
  MY_TYPE freq_resolution = (MY_TYPE)mydata->fs / n_fft;
  int f0_bin = (int)(f0 / freq_resolution + 0.5);  // Round to nearest bin

  // Question 16: Extract harmonics
  // Calculate magnitude and normalize (FFT output is scaled by n_fft for a cosine)
  int num_harmonics = 0;
  MY_TYPE* harmonic_freqs = new MY_TYPE[mydata->max_harmonics];
  MY_TYPE* harmonic_amps = new MY_TYPE[mydata->max_harmonics];
  MY_TYPE* harmonic_phases = new MY_TYPE[mydata->max_harmonics];

  for (int h = 1; h <= mydata->max_harmonics; h++) {
    int harmonic_bin = f0_bin * h;
    if (harmonic_bin < n_fft / 2) {
      harmonic_freqs[num_harmonics] = h * f0;

      // Calculate magnitude from real and imaginary parts
      MY_TYPE real_part = fft_real[harmonic_bin];
      MY_TYPE imag_part = fft_imag[harmonic_bin];
      MY_TYPE magnitude = sqrt(real_part * real_part + imag_part * imag_part);

      // Normalize: FFT gives DFT * n_fft for a real signal peak
      harmonic_amps[num_harmonics] = magnitude * 2.0 / n_fft;

      // Calculate phase
      harmonic_phases[num_harmonics] = atan2(imag_part, real_part);

      num_harmonics++;
    } else {
      break;
    }
  }

  // Question 17-18: Additive synthesis
  // Initialize output buffer to zero
  for (unsigned int i = 0; i < nBufferFrames; i++) {
    out[i] = 0.0;
  }

  // Synthesize each harmonic
  for (int h = 0; h < num_harmonics; h++) {
    MY_TYPE freq = harmonic_freqs[h];
    MY_TYPE amp = harmonic_amps[h];
    MY_TYPE phase = harmonic_phases[h];

    for (unsigned int i = 0; i < nBufferFrames; i++) {
      // Question 17: Simple synthesis (causes clicks due to phase discontinuity)
      // out[i] += amp * cos(2.0 * PI * freq * i / mydata->fs);

      // Question 18: Include phase for continuity between frames
      // Using sampleCounter to maintain global time reference
      // Phase = 2*pi*f*t + initial_phase (from FFT)
      // t = (sampleCounter + i) / fs
      MY_TYPE t = (MY_TYPE)(mydata->sampleCounter + i) / mydata->fs;
      MY_TYPE instantaneous_phase = 2.0 * PI * freq * t + phase;

      out[i] += amp * cos(instantaneous_phase);
    }
  }

  // Clamp output to prevent clipping
  for (unsigned int i = 0; i < nBufferFrames; i++) {
    if (out[i] > 1.0) out[i] = 1.0;
    if (out[i] < -1.0) out[i] = -1.0;
  }

  // Update sample counter for phase continuity
  mydata->sampleCounter += nBufferFrames;

  // Clean up
  delete[] fft_real;
  delete[] fft_imag;
  delete[] harmonic_freqs;
  delete[] harmonic_amps;
  delete[] harmonic_phases;

  // Question 11: Write output buffer to dump buffer
  write_buff_dump(out, nBufferFrames, mydata->dumpBufferOutput, mydata->dumpBufferSize, &mydata->dumpIndexOutput);

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

  unsigned int bufferFrames = 2048;
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

  // Question 15-18: Initialize FFT and synthesis parameters
  mydata->n_fft = get_nextpow2(bufferFrames);  // Next power of 2 >= bufferFrames
  mydata->max_harmonics = 50;  // Track up to 50 harmonics
  mydata->phases = new MY_TYPE[mydata->max_harmonics];
  for (int i = 0; i < mydata->max_harmonics; i++) {
    mydata->phases[i] = 0.0;
  }
  mydata->sampleCounter = 0;

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
  delete[] mydata->phases;
  delete mydata;

  return 0;
}
