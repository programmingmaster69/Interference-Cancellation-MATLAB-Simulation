# Interference-Cancellation-MATLAB-Simulation
* Introduction
  * With the evolution of wireless technologies such as 5G and Wi-Fi 6, MIMO-Orthogonal Frequency Division Multiplexing (MIMO-OFDM) systems have become popular due to their ability to support high data rates and robust signal transmission in multipath environments.
  * However, these systems are vulnerable to unintentional interference or jamming attacks, which can damage the communication by corrupting received signals.
Reactive Jamming is a threat to this system due to its adaptive nature—transmitting interference only when legitimate communication is detected. This makes it more energy-efficient and harder to detect than barrage jamming. To counter such attacks, modern receivers must be equipped with advanced interference cancellation (IC) and signal enhancement mechanisms.

* Objective
  * The objective of this project is to simulate a MIMO-OFDM system using MATLAB under jamming attacks and implement advanced signal processing techniques to detect, estimate, and mitigate the effects of jamming. The simulation is conducted using MATLAB, focusing on two key defenses:
    * Interference Cancellation (IC) using pilot-based channel estimation.
    * Sender Signal Enhancement (SSE) via optimal signal rotation in the transmitter's antenna space.

* System Design
We model a 2×2 MIMO-OFDM communication system. The OFDM frame includes:
  * 64 available subcarriers (48 used for data transmission),
  * 14 OFDM symbols per packet,
  * Cyclic prefix length - 16,
  * BPSK modulation/demodulation,
  * Additive White Gaussian Noise (AWGN),
  * Rayleigh fading for sender and jammer channels.
  * SNR (Signal to Noise Ratio): 16 dB

* To quantify system performance, we use:
  * Packet Delivery Rate (PDR): Percentage of packets correctly received (i.e., all symbols correctly demodulated).
  * Jam-to-Signal Ratio (JSR): Measured from the pilot residuals.
  * Comparison of three situations:
      * No countermeasure(defenseless),
      * With IC only,
      * With IC + SSE.

* Results:
  PDR plots were generated for each configuration across increasing jamming power levels, clearly showing the gain from IC and SSE.
  Jamming power was varied between 0 - 10 W.
  






<img width="1000" height="479" alt="image" src="https://github.com/user-attachments/assets/1932ad87-9eb0-4413-8ea7-989103bc8f55" />
