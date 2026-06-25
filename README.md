# Water-Level-Control-of-Three-Interconnected-Tanks

This project models and controls a three‑tank water level system using nonlinear dynamics, linearization, and both linear and nonlinear Model Predictive Control (MPC). The workflow includes deriving the nonlinear ODE model, linearizing it around an equilibrium point and designing four MPC controllers (linear and nonlinear, with and without embedded constraints), each tested under nominal and disturbed conditions.

## System Modeling
- Deriving three coupled nonlinear ODEs describing water levels and inter‑tank flows including system disturbances
- Calculation of the equilibrium point based on the tank water level setpoints
- System linearization around the operating-equilibrium point
- Using discrete‑time state‑space representation for MPCs implementation

## Controllers Implemented
- Linear MPC based on the linearized model, without embedded constraints (constraints handled externally), tested with and without disturbances
- Linear MPC based on the linearized model, with embedded constraints (constraints handled externally), tested with and without disturbances
- Nonlinear MPC based on the nonlinear model, without embedded constraints (constraints handled externally), tested with and without disturbances
- Nonlinear MPC based on the nonlinear model, with embedded constraints (constraints handled externally), tested with and without disturbances

## Repository Contents
- `Linear_MPC_No_Constraints/`
  - `Linear_MPC_No_Constraints_No_Disturbances.m`
  - `Linear_MPC_No_Constraints_Disturbances.m`
- `Linear_MPC_Constraints/`
  - `Linear_MPC_Constraints_No_Disturbances.m`
  - `Linear_MPC_Constraints_Disturbances.m`
- `Non_Linear_MPC_No_Constraints/`
  - `Non_Linear_MPC_No_Constraints_No_Disturbances.m`
  - `Non_Linear_MPC_No_Constraints_Disturbances.m`
-`Non_Linear_MPC_Constraints/`
  - `Non_Linear_MPC_Constraints_No_Disturbances.m`
  - `Non_Linear_MPC_Constraints_Disturbances.m`
- `functions/`
  - `myStateFunction.m` - Nonlinear state-space model
  - `myStateJacobian.m` Jacobian of nonlinear ODEs with respect to states and inputs
  - `obj_fun.m` - Objective function for the linear MPC with embedded constraints
  - `ode_fun.m` - Nonlinear model

## How to Run
- Choose the controller script
- Run it directly in MATLAB
- The script automatically:
  - loads the model and parameters
  - runs the MPC simulation
  - applies constraints (embedded or external)
  - generates plots of tank levels and control inputs
