# ESRA Motor Performance Script
This project was used during my senior year of college to help calculate key performance parameters of a solid rocket motor. It calculates maximum pressure, and maximum thrust. In addition, the program calculates delivered and theoretical characteristic velocity, impulse, specific impulse, nozzle exit pressure, and coefficient of thrust.
# Input Constraints
* CSV file contains three columns in the following order: time, pressure, load
* CSV file must be titled with format "~X.XXXXX\~X.XXXXX\~X.XXXXX.csv" where each corresponding value represents: "Nozzle throat diameter [in.] ~ Nozzle exit diameter [in.] ~ mass of propellant [lbs.]".
* Theoretical specific heat ratio, combustion temperature, and molecular weight are specific to LS propellant. 
* Calibration slopes and intercepts
