# Coherence-Based Source Counter
*Counts the number of speech sources during a meeting using two microphone nodes*

[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.879279-blue.svg?style=flat-square)](https://doi.org/10.5281/zenodo.879279)
[![GitHub release](https://img.shields.io/github/release/JacobD10/CoherenceBasedSourceCounter/all.svg?style=flat-square)](https://github.com/JacobD10/CoherenceBasedSourceCounter/releases)
[![GitHub commits](https://img.shields.io/github/commits-since/JacobD10/CoherenceBasedSourceCounter/1.0.0.svg?style=flat-square)](https://github.com/JacobD10/CoherenceBasedSourceCounter/commits/master)
[![Hit Count](https://hitt.herokuapp.com/JacobD10/CoherenceBasedSourceCounter.svg?style=flat-square)](https://github.com/JacobD10/CoherenceBasedSourceCounter)
[![Github file size](https://reposs.herokuapp.com/?path=JacobD10/CoherenceBasedSourceCounter&style=flat-square)](https://github.com/JacobD10/CoherenceBasedSourceCounter/archive/master.zip)
[![license](https://img.shields.io/github/license/JacobD10/CoherenceBasedSourceCounter.svg?style=flat-square)](https://github.com/JacobD10/CoherenceBasedSourceCounter/blob/master/LICENSE)
[![Twitter URL](https://img.shields.io/twitter/url/http/shields.io.svg?style=social)](https://twitter.com/intent/tweet?url=https%3A%2F%2Fgithub.com%2FJacobD10%2FCoherenceBasedSourceCounter&via=_JacobDonley&text=Check%20out%20the%20Two-Microphone%20Source%20Counter%20for%20%23MATLAB%21&hashtags=software%20%23code%20%23audio)

## **Reference**
S. Pasha, J. Donley and C. Ritz, "**Blind Speaker Counting in Highly Reverberant Environments by Clustering Coherence Features,**" in *Asia-Pacific Signal & Information Processing Association Annual Summit and Conference (APSIPA ASC)*.  IEEE, Dec. 2017, pp.1-5.

---

## Description
This MATLAB function clusters the Magnitude-Squared Coherence (MSC) of two channel signals and uses k-means to estimate the number of speech sources contributing to the input signal.

## First Use and Example
**First Run**: In order to generate test signals for the source counting algorithm, a room impulse response (RIR) is generated from the [Room Impulse Response Generator by Emanuel A.P. Habets](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator). Please download, compile and ensure that the function called `rir_generator()` is on the MATLAB search path prior to running `example.m` from this repository.

1. [Download the latest release of this repository](https://github.com/JacobD10/CoherenceBasedSourceCounter/releases) and navigate to the directory in MATLAB.
2. Open [`example.m`](https://github.com/JacobD10/CoherenceBasedSourceCounter/blob/master/example.m) and set `SynthParams.SpeechDir` to the directory which containes all the speech (talker) audio files to evaluate with.
Note: If you have a copy of the [TIMIT corpus](https://catalog.ldc.upenn.edu/ldc93s1) you may use the helper function [`ConcatTIMITtalkers()`](https://github.com/JacobD10/CoherenceBasedSourceCounter/blob/master/ConcatTIMITtalkers.m) to generate the talker audio files.
3. Run `example.m`.

## Function Usage
The function which implements the source counting method is located in [`cbsc.m`](https://github.com/JacobD10/CoherenceBasedSourceCounter/blob/master/cbsc.m). Help and an example of how to use the function are located in the header of the file. The function requires only two microphone signals when using the default `AnalysisParams`. The microphone signals should be of the size `[samples microphones nodes]` (_row_, _column_, _page_) where there are only two columns (two microphones).

## Contact 
Shahab Pasha - sp900@uowmail.edu.au

Jacob Donley - jrd089@uowmail.edu.au

Christian Ritz - critz@uow.edu.au

