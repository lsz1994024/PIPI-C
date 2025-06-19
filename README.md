# PIPI-C
Identifying PTM Combinations in Peptides Using Mixed Integer Linear Programming

Contact: slaiad@connect.ust.hk
https://bioinformatics.hkust.edu.hk/
 
## Requirement
- Java 1.8.
Download and install JDK 8 from https://www.oracle.com/java/technologies/downloads/#java8
- Gurobi 10.0.3.
Download and install Gurobi from https://www.gurobi.com/downloads/

Usage:

with Java 8 (recommended):
java -Xmx8g -jar PIPI.jar parameter.def

with Java of higher versions:
java -Xmx8g -jar --add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED PIPI.jar parameter.def




