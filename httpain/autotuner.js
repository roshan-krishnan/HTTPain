/* AUTOTUNER.JS */

// initial variables
var M_PI = 3.14159265358979323846;              // just a constant for pi
var halfStep = 17.32 / 16.35;                   // multiplying value to shift a frequency up by a half-step
var bufferSize = 256;                           // the buffer size for our script processor node
var osamp = 16;                                 // oversampling factor
var fftSize = bufferSize*osamp;                 // our FFT frame size
var fftSize2 = fftSize/2;                       // half the FFT size, useful for spectral processing
var hopSize = bufferSize;                       // the hop size for our overlap add- just use SPN buffer size
var fs;                                         // sample rate
var fifoInCounter = 0;                          // counter to step through the constant FIFO input buffer
var fifoOutCounter = 0;                         // counter to step through the constant FIFO output buffer
var expPhaseDiff = 2.0*M_PI*hopSize/fftSize;    // the expected phase difference between 2 FFT bins
var freqPerBin;                                 // frequency per bin. will just be fs / FFT size
var rec;                                        // variable to hold a RecorderJS object.
var isMaj = true;                               // flag that denotes whether we've selected a major key
var keyID;                                      // ID of a key- C is 1, C# is 2, ...., B is 12.
var isRecording = false;                        // flag that indicates whether a recording is in progress

// Basic note frequency lookup table for a MAJOR scale from C0 to B8.
var majorScale = [
    16.35,18.35,20.60,21.83,24.50,27.50,30.87,32.70,
    36.71,41.20,43.65,49.00,55.00,61.74,65.41,73.42,
    82.41,87.31,98.00,110.00,123.47,130.81,146.83,164.81,
    174.61,196.00,220.00,246.94,261.63,293.66,329.63,349.23,
    392.00,440.00,493.88,523.25,587.33,659.25,698.46,783.99,
    880.00,987.77,1046.50,1174.66,1318.51,1396.91,1567.98,1760.00,
    1975.53,2093.00,2349.32,2637.02,2793.83,3135.96,3520.00,3951.07,
    4186.01,4698.63,5274.04,5587.65,6271.93,7040.00,7902.13
];

// Basic note frequency lookup table for a NATURAL MINOR scale from C0 to B8.
var minorScale = [
    16.35,18.35,19.45,21.83,24.50,25.96,29.14,32.70,
    36.71,38.89,43.65,49.00,51.91,58.27,65.41,73.42,
    77.78,87.31,98.00,103.83,116.54,130.81,146.83,155.56,
    174.61,196.00,207.65,233.08,261.63,293.66,311.13,349.23,
    392.00,415.30,466.16,523.25,587.33,622.25,698.46,783.99,
    830.61,932.33,1046.50,1174.66,1244.51,1396.91,1567.98,1661.22,
    1864.66,2093.00,2349.32,2489.02,2793.83,3135.96,3322.44,3729.31,
    4186.01,4698.63,4978.03,5587.65,6271.93,6644.88,7458.62
];

// array of pictures of T-Pain, because why not?
var tPainArray = [
    "../images/tpain1.jpg","../images/tpain2.png","../images/tpain3.jpg","../images/tpain4.png","../images/tpain5.jpg",
    "../images/tpain6.jpg","../images/tpain7.jpg","../images/tpain8.jpg","../images/tpain9.jpg","../images/tpain10.png"
];


// wait for the window to load
window.onload = function() {
    // get the audio context
    window.AudioContext = window.AudioContext || window.webkitAudioContext;
    var context = new AudioContext();
    fs = context.sampleRate;
    freqPerBin = fs/fftSize;
    
    // select a random picture of T-Pain when
    randomTPain();
    
    // set the text indicating our current key (initializes to C Major)
    keyText(keyID, isMaj);
    
    // create an FFT object
    var fft = new FFT(fftSize, fs);
    
    // create our input/output windowing buffers
    var fifoIn = new Float32Array(fftSize);
    var fifoOut = new Float32Array(fftSize);
    
    // create static arrays / vars
    var aMagnitude = new Float32Array(fftSize);
    var aFrequency = new Float32Array(fftSize);
    var lastPhase = new Float32Array(fftSize2+1);
    var sumPhase = new Float32Array(fftSize2+1);
    
    // make sure we have getUserMedia
    navigator.getUserMedia = navigator.getUserMedia ||
                         navigator.webkitGetUserMedia ||
                         navigator.mozGetUserMedia;

    if(navigator.getUserMedia){
        // start streaming from microphone
        navigator.getUserMedia({audio: true}, function(stream) {
            var microphone = context.createMediaStreamSource(stream);
            var processor = context.createScriptProcessor(bufferSize, 1, 1);
            
            // there's not much vocal energy above 8k, so filter it out
            var lpfIn = context.createBiquadFilter();
            lpfIn.type = "lowpass";
            lpfIn.frequency.value = 8000;
            var lpfOut = context.createBiquadFilter();
            lpfOut.type = "lowpass";
            lpfOut.frequency.value = 8000;
            
            // wire our nodes up
            microphone.connect(lpfIn);
            lpfIn.connect(processor);
            processor.connect(lpfOut);
            lpfOut.connect(context.destination);
            
            // create a new RecorderJS object
            rec = new Recorder(processor, {workerPath: "httpain/recorderWorker.js"});
            
            // do buffered audio processing
            processor.onaudioprocess = function(audioData){
                // grab the input and output buffers
                var inputData = audioData.inputBuffer.getChannelData(0);
                var outputData = audioData.outputBuffer.getChannelData(0);
                
                // write the oldest frame in the fifo out buffer to output
                // and clear that space in the fifo out buffer
                for(var i = 0; i < hopSize; i++){
                    outputData[i] = fifoOut[fifoInCounter+i];
                    fifoOut[fifoInCounter+i] = 0.0;
                }
                
                // write the current frame into my input buffer
                for(var i = 0; i < bufferSize; i++)
                    fifoIn[fifoInCounter+i] = inputData[i];

                // make sure we do windowing in order
                var next = fifoInCounter + hopSize;
                if(next == fftSize)
                    next = 0;
                
                // do the actual windowing
                var windowed = hannWindow(fifoIn, next);

                // forward FFT
                fft.forward(windowed);
                var real = fft.real;
                var imag = fft.imag;
                var phase = 0;
                var temp = 0;
                
                // Fourier analysis
                for(var i = 0; i <= fftSize2; i++){
                    // find the magnitude and phase
                    aMagnitude[i] = 2.0*Math.sqrt(real[i]*real[i] + imag[i]*imag[i]);
                    phase = Math.atan2(imag[i], real[i]);
                    
                    // find the phase difference
                    temp = phase - lastPhase[i];
                    lastPhase[i] = phase;
                    
                    // subtract the expected difference
                    temp -= i*expPhaseDiff;
                    
                    // map the phase difference between +/- pi
                    var qpd = temp/M_PI;
                    if(qpd >= 0)
                        qpd += (qpd & 1);
                    else
                        qpd -= (qpd & 1);
                    
                    // get deviation and compute the true frequency
                    temp -= M_PI*qpd;
                    temp = (osamp*temp)/(2.0*M_PI);
                    temp = i*freqPerBin + temp*freqPerBin;
                    aFrequency[i] = temp;  
                }

                // for now we just get the frequency with the highest magnitude.
                // that should be good enough for voice, but it should be noted
                // that the the frequency with the highest magnitude isn't always
                // the fundamental, and in the future it would be wise to implement
                // a pitch detection algorithm here.
                var maxFreqIndex = 0;
                for(var i = 0; i < aMagnitude.length; i++){
                    if(aMagnitude[i] > maxFreqIndex)
                        maxFreqIndex = i;
                }

                // get the lookup table of note frequencies for the currently selected key.
                var noteArray = multKey(keyID, isMaj);
                
                // get the shift factor
                var shiftFactor = getShiftFactor(aFrequency[maxFreqIndex], noteArray);
                
                // create synthesis arrays
                var sMagnitude = new Float32Array(fftSize);
                var sFrequency = new Float32Array(fftSize);
                var index;
                
                // do the actual pitch shifting
                for(var i = 0; i <= fftSize2; i++){
                    index = Math.round(i*shiftFactor);
                    if(index <= fftSize2){
                        sMagnitude[index] += aMagnitude[i];
                        sFrequency[index] = aFrequency[i]*shiftFactor;
                    }
                }
                
                // create a temp variable to hold magnitude in the synthesis loop
                var magn;
                
                // do synthesis
                for(var i = 0; i <= fftSize2; i++){
                    // get magnitude and freq values
                    magn = sMagnitude[i];
                    temp = sFrequency[i];
                    
                    // reverse our frequency bin deviation calculations
                    temp -= i*freqPerBin;
                    temp /= freqPerBin;
                    temp = 2.0*M_PI*temp/osamp;
                    
                    // factor the phase back in
                    temp += i*expPhaseDiff;
                    sumPhase[i] += temp;
                    phase = sumPhase[i];
                    
                    // store it away
                    real[i] = magn*Math.cos(phase);
                    imag[i] = magn*Math.sin(phase);
                }
                
                // zero negative frequencies. we only want the real result
                for(var i = fftSize2+1; i < fftSize; i++){
                    real[i] = 0.0;
                    imag[i] = 0.0;
                }

                // inverse FFT
                var ifft = fft.inverse(real, imag);

                // windowing
                var processed = hannWindow(ifft, 0);
                
                // add result to fifo out in time order
                var inc = fifoOutCounter;
                for(var i = 0; i < fftSize; i++){
                    fifoOut[i] += processed[inc]/osamp;
                    inc++;
                    if(inc == fftSize)
                        inc = 0;
                }

                // decrement the fifo out frame counter and wrap if necessary
                fifoOutCounter -= hopSize;
                if(fifoOutCounter == -hopSize)
                    fifoOutCounter = fftSize - hopSize;
                
                // increment the fifo in frame counter and wrap if necessary
                fifoInCounter += hopSize;
                if(fifoInCounter == fftSize)
                    fifoInCounter = 0;

            }
        },
        function(err){
            console.log('Streaming error!');
        });
    }else{
        console.log('getUserMedia not supported!');
    }
}

// takes an array of data and an index (COUNTER) at which to start windowing, and implements
// a Hann window. Note that DATA is maintained and we return a brand new array.
function hannWindow(data, counter){
    var dataCopy = new Float32Array(data.length);
    for(var i = 0; i < dataCopy.length; i++){
        var mult = 0.5*(1.0 - Math.cos(2.0*M_PI*i/data.length));
        dataCopy[i] = data[counter]*mult;
        counter++;
        if(counter >= dataCopy.length)
            counter = 0;
    }
    return dataCopy;
}

// takes in our input frequency VALUE and an array NOTES containing all the note frequencies
// in the desired key, and returns the closest value in NOTES to VALUE.
function getClosestVal(value, notes){
    var curr = notes[0];
    for(var i = 0; i < notes.length; i++){
        var c1 = Math.abs(value - notes[i]);
        var c2 = Math.abs(value - curr);
        if(c1 < c2)
            curr = notes[i];
    }
    return curr;
}

// calls the getClosestValue helper function to get the closest value to FREQVALUE,
// then computes and returns the ratio between the two.
function getShiftFactor(freqValue, notes){
    var closest = getClosestVal(freqValue, notes);
    var ratio = closest/freqValue;
    return ratio;
}

// starts recording if a recording is not already in progress.
function startRecord(){
    if(!isRecording){
        rec.record();
        isRecording = true;
    }
}

// stops recording if a recording is in progress, and exports and downloads the WAV.
function stopRecord(){
    if(isRecording){
        rec.stop();
        rec.exportWAV(function(blob){
            Recorder.forceDownload(blob);
        });
        rec.clear();
        isRecording = false;
    }
}

// multiplies the major or minor (depending on the value of boolean ISMAJ) scale tables
// by VALUE, then returns them in a new array.
function multByVal(value, isMaj){
    var outputArr = [];
    if(isMaj){
        for(var i = 0; i < majorScale.length; i++)
            outputArr.push(majorScale[i]*value);
    }else{
        for(var i = 0; i < minorScale.length; i++)
            outputArr.push(minorScale[i]*value);
    }
    return outputArr;
}

// uses MULTBYVAL(VALUE, ISMAJ) to multiply the template scales according to the musical key
// signified by ID and the major/minor quality indicated by ISMAJ. Returns an array of the notes
// in the desired scale.
function multKey(id, isMaj){
    var notes;
    switch(id){
        case "1":
            notes = multByVal(1, isMaj);
            break;
        case "2":
            notes = multByVal(halfStep, isMaj);
            break;
        case "3":
            notes = multByVal(Math.pow(halfStep, 2), isMaj);
            break;
        case "4":
            notes = multByVal(Math.pow(halfStep, 3), isMaj);
            break;
        case "5":
            notes = multByVal(Math.pow(halfStep, 4), isMaj);
            break;
        case "6":
            notes = multByVal(Math.pow(halfStep, 5), isMaj);
            break;
        case "7":
            notes = multByVal(Math.pow(halfStep, 6), isMaj);
            break;
        case "8":
            notes = multByVal(Math.pow(halfStep, 7), isMaj);
            break;
        case "9":
            notes = multByVal(Math.pow(halfStep, 8), isMaj);
            break;
        case "10":
            notes = multByVal(Math.pow(halfStep, 9), isMaj);
            break;
        case "11":
            notes = multByVal(Math.pow(halfStep, 10), isMaj);
            break;
        case "12":
            notes = multByVal(Math.pow(halfStep, 11), isMaj);
            break;
        default:
            notes = multByVal(1, isMaj);
            break;
    }
    return notes;
}

// sets the key based on which button is pressed.
function setKey(elem){
    keyID = elem.id;
    randomTPain();
    keyText(keyID, isMaj);
}

// sets major or minor based on which button is pressed.
function majMin(elem){
    if(elem.id == "maj"){
        isMaj = true;
    }else{
        isMaj = false;
    }
    randomTPain();
    keyText(keyID, isMaj);
}

// generates a random picture of T-Pain from our array of images.
function randomTPain(){
    var randIndex = Math.floor(Math.random() * 10);
    var text = "<img src=" + tPainArray[randIndex] + "></img>";
    document.getElementById("tpain").innerHTML = text;
}

// generates text that displays the current key based on the values
// of ID (key) and ISMAJ (major/minor quality).
function keyText(id, isMaj){
    var str;
    switch(id){
        case "1":
            str = "C";
            break;
        case "2":
            str = "C#/Db";
            break;
        case "3":
            str = "D";
            break;
        case "4":
            str = "D#/Eb";
            break;
        case "5":
            str = "E";
            break;
        case "6":
            str = "F";
            break;
        case "7":
            str = "F#/Gb";
            break;
        case "8":
            str = "G";
            break;
        case "9":
            str = "G#/Ab";
            break;
        case "10":
            str = "A";
            break;
        case "11":
            str = "A#/Bb";
            break;
        case "12":
            str = "B";
            break;
        default:
            str = "C";
            break;
    }
    if(isMaj)
        str += " Major";
    else
        str += " Minor";
    document.getElementById("keytext").innerHTML = "Key: "+str;
}

// Copyright (c) 2010 Corban Brook, released under the MIT license
// Fourier Transform Module used by DFT, FFT, RFT
// Slightly modified for packed DC/nyquist...

function FourierTransform(bufferSize, sampleRate) {
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;
  this.bandwidth  = 2 / bufferSize * sampleRate / 2;

  this.spectrum   = new Float32Array(bufferSize/2);
  this.real       = new Float32Array(bufferSize);
  this.imag       = new Float32Array(bufferSize);

  this.peakBand   = 0;
  this.peak       = 0;

  /**
   * Calculates the *middle* frequency of an FFT band.
   *
   * @param {Number} index The index of the FFT band.
   *
   * @returns The middle frequency in Hz.
   */
  this.getBandFrequency = function(index) {
    return this.bandwidth * index + this.bandwidth / 2;
  };

  this.calculateSpectrum = function() {
    var spectrum  = this.spectrum,
        real      = this.real,
        imag      = this.imag,
        bSi       = 2 / this.bufferSize,
        sqrt      = Math.sqrt,
        rval, 
        ival,
        mag;

    for (var i = 0, N = bufferSize/2; i < N; i++) {
      rval = real[i];
      ival = imag[i];
      mag = bSi * sqrt(rval * rval + ival * ival);

      if (mag > this.peak) {
        this.peakBand = i;
        this.peak = mag;
      }

      spectrum[i] = mag;
    }
  };
}

/**
 * FFT is a class for calculating the Discrete Fourier Transform of a signal
 * with the Fast Fourier Transform algorithm.
 *
 * @param {Number} bufferSize The size of the sample buffer to be computed. Must be power of 2
 * @param {Number} sampleRate The sampleRate of the buffer (eg. 44100)
 *
 * @constructor
 */
function FFT(bufferSize, sampleRate) {
  FourierTransform.call(this, bufferSize, sampleRate);
   
  this.reverseTable = new Uint32Array(bufferSize);

  var limit = 1;
  var bit = bufferSize >> 1;

  var i;

  while (limit < bufferSize) {
    for (i = 0; i < limit; i++) {
      this.reverseTable[i + limit] = this.reverseTable[i] + bit;
    }

    limit = limit << 1;
    bit = bit >> 1;
  }

  this.sinTable = new Float32Array(bufferSize);
  this.cosTable = new Float32Array(bufferSize);

  for (i = 0; i < bufferSize; i++) {
    this.sinTable[i] = Math.sin(-Math.PI/i);
    this.cosTable[i] = Math.cos(-Math.PI/i);
  }
}

/**
 * Performs a forward tranform on the sample buffer.
 * Converts a time domain signal to frequency domain spectra.
 *
 * @param {Array} buffer The sample buffer. Buffer Length must be power of 2
 *
 * @returns The frequency spectrum array
 */
FFT.prototype.forward = function(buffer) {
  // Locally scope variables for speed up
  var bufferSize      = this.bufferSize,
      cosTable        = this.cosTable,
      sinTable        = this.sinTable,
      reverseTable    = this.reverseTable,
      real            = this.real,
      imag            = this.imag,
      spectrum        = this.spectrum;

  var k = Math.floor(Math.log(bufferSize) / Math.LN2);

  if (Math.pow(2, k) !== bufferSize) { throw "Invalid buffer size, must be a power of 2."; }
  if (bufferSize !== buffer.length)  { throw "Supplied buffer is not the same size as defined FFT. FFT Size: " + bufferSize + " Buffer Size: " + buffer.length; }

  var halfSize = 1,
      phaseShiftStepReal,
      phaseShiftStepImag,
      currentPhaseShiftReal,
      currentPhaseShiftImag,
      off,
      tr,
      ti,
      tmpReal,
      i;

  for (i = 0; i < bufferSize; i++) {
    real[i] = buffer[reverseTable[i]];
    imag[i] = 0;
  }

  while (halfSize < bufferSize) {
    //phaseShiftStepReal = Math.cos(-Math.PI/halfSize);
    //phaseShiftStepImag = Math.sin(-Math.PI/halfSize);
    phaseShiftStepReal = cosTable[halfSize];
    phaseShiftStepImag = sinTable[halfSize];
    
    currentPhaseShiftReal = 1;
    currentPhaseShiftImag = 0;

    for (var fftStep = 0; fftStep < halfSize; fftStep++) {
      i = fftStep;

      while (i < bufferSize) {
        off = i + halfSize;
        tr = (currentPhaseShiftReal * real[off]) - (currentPhaseShiftImag * imag[off]);
        ti = (currentPhaseShiftReal * imag[off]) + (currentPhaseShiftImag * real[off]);

        real[off] = real[i] - tr;
        imag[off] = imag[i] - ti;
        real[i] += tr;
        imag[i] += ti;

        i += halfSize << 1;
      }

      tmpReal = currentPhaseShiftReal;
      currentPhaseShiftReal = (tmpReal * phaseShiftStepReal) - (currentPhaseShiftImag * phaseShiftStepImag);
      currentPhaseShiftImag = (tmpReal * phaseShiftStepImag) + (currentPhaseShiftImag * phaseShiftStepReal);
    }

    halfSize = halfSize << 1;
  }
  
  // Pack nyquist component.
  imag[0] = real[bufferSize / 2];
};

FFT.prototype.inverse = function(real, imag) {
  // Locally scope variables for speed up
  var bufferSize      = this.bufferSize,
      cosTable        = this.cosTable,
      sinTable        = this.sinTable,
      reverseTable    = this.reverseTable,
      spectrum        = this.spectrum;
     
      real = real || this.real;
      imag = imag || this.imag;

  var halfSize = 1,
      phaseShiftStepReal,
      phaseShiftStepImag,
      currentPhaseShiftReal,
      currentPhaseShiftImag,
      off,
      tr,
      ti,
      tmpReal,
      i;

      // Unpack and create mirror image.
      // This isn't that efficient, but let's us avoid having to deal with the mirror image part
      // when processing.
      var n = bufferSize;
      var nyquist = imag[0];
      imag[0] = 0;
      real[n / 2] = nyquist;
      imag[n / 2] = 0;

      // Mirror image complex conjugate.
      for (i = 1 + n / 2; i < n; i++) {
          real[i] = real[n - i];
          imag[i] = -imag[n - i];
      }

  for (i = 0; i < bufferSize; i++) {
    imag[i] *= -1;
  }

  var revReal = new Float32Array(bufferSize);
  var revImag = new Float32Array(bufferSize);

  
  
  
   
  for (i = 0; i < real.length; i++) {
    revReal[i] = real[reverseTable[i]];
    revImag[i] = imag[reverseTable[i]];
  }
   
  real = revReal;
  imag = revImag;
  
  while (halfSize < bufferSize) {
    phaseShiftStepReal = cosTable[halfSize];
    phaseShiftStepImag = sinTable[halfSize];
    currentPhaseShiftReal = 1;
    currentPhaseShiftImag = 0;
  
    for (var fftStep = 0; fftStep < halfSize; fftStep++) {
      i = fftStep;
  
      while (i < bufferSize) {
        off = i + halfSize;
        tr = (currentPhaseShiftReal * real[off]) - (currentPhaseShiftImag * imag[off]);
        ti = (currentPhaseShiftReal * imag[off]) + (currentPhaseShiftImag * real[off]);
  
        real[off] = real[i] - tr;
        imag[off] = imag[i] - ti;
        real[i] += tr;
        imag[i] += ti;
  
        i += halfSize << 1;
      }
  
      tmpReal = currentPhaseShiftReal;
      currentPhaseShiftReal = (tmpReal * phaseShiftStepReal) - (currentPhaseShiftImag * phaseShiftStepImag);
      currentPhaseShiftImag = (tmpReal * phaseShiftStepImag) + (currentPhaseShiftImag * phaseShiftStepReal);
    }
  
    halfSize = halfSize << 1;
  }
  
  var buffer = new Float32Array(bufferSize); // this should be reused instead
  for (i = 0; i < bufferSize; i++) {
    buffer[i] = real[i] / bufferSize;
  }

  return buffer;
};

(function(window){

  var WORKER_PATH = 'recorderWorker.js';

  var Recorder = function(source, cfg){
    var config = cfg || {};
    var bufferLen = config.bufferLen || 4096;
    var numChannels = config.numChannels || 2;
    this.context = source.context;
    this.node = (this.context.createScriptProcessor ||
                 this.context.createJavaScriptNode).call(this.context,
                 bufferLen, numChannels, numChannels);
    var worker = new Worker(config.workerPath || WORKER_PATH);
    worker.postMessage({
      command: 'init',
      config: {
        sampleRate: this.context.sampleRate,
        numChannels: numChannels
      }
    });
    var recording = false,
      currCallback;

    this.node.onaudioprocess = function(e){
      if (!recording) return;
      var buffer = [];
      for (var channel = 0; channel < numChannels; channel++){
          buffer.push(e.inputBuffer.getChannelData(channel));
      }
      worker.postMessage({
        command: 'record',
        buffer: buffer
      });
    }

    this.configure = function(cfg){
      for (var prop in cfg){
        if (cfg.hasOwnProperty(prop)){
          config[prop] = cfg[prop];
        }
      }
    }

    this.record = function(){
      recording = true;
    }

    this.stop = function(){
      recording = false;
    }

    this.clear = function(){
      worker.postMessage({ command: 'clear' });
    }

    this.getBuffer = function(cb) {
      currCallback = cb || config.callback;
      worker.postMessage({ command: 'getBuffer' })
    }

    this.exportWAV = function(cb, type){
      currCallback = cb || config.callback;
      type = type || config.type || 'audio/wav';
      if (!currCallback) throw new Error('Callback not set');
      worker.postMessage({
        command: 'exportWAV',
        type: type
      });
    }

    worker.onmessage = function(e){
      var blob = e.data;
      currCallback(blob);
    }

    source.connect(this.node);
    this.node.connect(this.context.destination);    //this should not be necessary
  };

  Recorder.forceDownload = function(blob, filename){
    var url = (window.URL || window.webkitURL).createObjectURL(blob);
    var link = window.document.createElement('a');
    link.href = url;
    link.download = filename || 'output.wav';
    var click = document.createEvent("Event");
    click.initEvent("click", true, true);
    link.dispatchEvent(click);
  }

  window.Recorder = Recorder;

})(window);