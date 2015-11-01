# HTTPain
AutoTune for the browser, feat. everyone's favorite AutoTune popularizer.

I used <a href=https://github.com/corbanbrook/dsp.js/>DSP.js</a> for my FFT/IFFT,
<a href=https://github.com/mattdiamond/Recorderjs>Recorder.js</a> for the recording/WAV export functionality, 
and Stephan Bernsee's excellent article <a href=http://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/>Pitch Shifting 
Using the Fourier Transform</a> as inspiration.

There's still some work to be done. A little bit of aliasing funny stuff in the higher frequencies, and it would be better to use a real pitch-detection algorithm (maybe HPS?) rather than just assuming the bin with the largest magnitude contains the fundamental. But hey, if you've ever wanted to sing "Buy U A Drank" alone in your room at your computer, well, now you can.
