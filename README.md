# DSB_Modulation

This is a DSB-LC (Double Sideband Large Carrier) Modulation model.

Modulator Description:

-----------------------Originatopr Side---------------------------			Free Space				-------------------Recipient Side---------------------------
------------------------------------------------------------------									------------------------------------------------------------
-------v_m----->  * v_c ------> * kAM ---------> + v_c ------------> ---------v_mod---------------> conv with BPF ------> Envelope Detector---------------------
------------------------------------------------------------------									------------------------------------------------------------


v_m - Data Signal (m for message)

v_c - Carrier Signal (cosine with very high frequency)

v_mod - modulated Signal

kAM - modulation index

Model Flow:

1.	The modulated signal is generated from the carrier signal and the mesage signal by the following system:
    * 	v_mod = (1 + kAM*v_m) * v_c = v_c + kAM*v_c*v_m
	This transformation is adding a DC to the signal and shifting it to very high frequency.
	Why is it serving us?
	DC Addition - so the recipient side will be able to use an envelope detector (takes v_c as an input).
	Shifting - many channels can operate in the same space and use different center frequency.
	

2. 	The modulated signal is passing the "free space" channel between the Originator and Recipient sides.
	As the signal passes the channel a gaussian noise being addded to the signal to reflect a real "free space" channel.
	In order to reflect the difference between reconstructed signals according to noise variance there were 3 different noises constructed.
	
3.  The Recipient side gets the signal after addition of the noise and filters it with Band Pass Filter with center frequency of the carrier.
	Than we use an envelope detector in order to extract the message signal from the filtered signal.

Repository files:

There are 3 files in this repository:

- "project1.m" 					- 	A matlab code file that needs to be run in a matlab enviroment and has no dependencies (you can run it as a single file).
- "DSB-LC Channel.png" 			- 	An image depicses the structure of a simple DSB-LC Modulator.
- "Frequency Description.png" 	- 	multiple graphs shows Modulation process in the frequency domain.
- "readme.md"					- 	Full Repository Description.

