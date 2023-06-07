# Integrated Lateral Stability and Traction Control for Electric Vehicles: Model Predictive Controller

## Overview

This GitHub repository contains the implementation of a Model Predictive Controller (MPC) designed for integrated lateral stability, traction/braking control, and rollover prevention of electric vehicles in very high-speed (VHS) racing applications. The controller aims to achieve stable control at speeds significantly greater than typical highway speed limits.

This controller has been developed as part of a research study that addresses the limitations of traditional non-MPC controllers used in high-speed racing scenarios. The paper associated with this repository presents the advantages of a state-of-the-art dynamic model that incorporates rollover prevention and linearizes the tire model to optimize computation time.

## Research Highlights

The key highlights of the research and the associated paper are as follows:

- Identification of the advantages of a dynamic model that integrates rollover prevention and linearized tire model in the MPC framework.
- Proposal of a novel model predictive controller for lateral stability control, specifically targeting stable control at top speeds significantly greater than typical highway speed limits.
- Testing and evaluation of the proposed controller in simulation environments associated with the Indy Autonomous Challenge, which includes real-world racing conditions such as road banking angles, lateral position tracking, and a different suspension model of the Dallara Indy Lights chassis.
- Promising results achieved with a low solver time in Python, reaching as low as 50 Hz, and a lateral error of 30 cm at speeds of 45 m/s.

## Usage

To utilize the Model Predictive Controller (MPC) implementation provided in this repository, please follow these steps:

1. Clone the repository to your local machine using the following command:
`git clone git@github.com:jadyahya/Roll-Yaw-and-Lateral-Velocity-MPC.git
`
2. Modify and configure the controller parameters to suit your specific vehicle and racing conditions.

3. Run the controller and evaluate its performance using the provided simulation environments and scenarios.


# Note: The file Class_CasADi_Implemented.py is run through the LGSVL simulator (open source) and the remainder of the software stack (localization, pathplanning, etc.). The latter are proprietary and not available in this repository.

## Contribution and Feedback

We believe that the implementation and findings presented in this repository can greatly benefit autonomous vehicle applicationa generally, and the racing community specifically, particularly those participating in high-speed racing events such as the Indy Autonomous Challenge. We encourage researchers, developers, and enthusiasts to contribute to this project and provide feedback to further improve the performance and capabilities of the Model Predictive Controller (MPC).

## Contact

For any questions, suggestions, or inquiries related to this repository, please contact the project maintainers at [jadraafatyahya@gmail.com].

We appreciate your interest in this research and look forward to collaborating with you to enhance the performance and safety of electric vehicles in high-speed racing applications.
