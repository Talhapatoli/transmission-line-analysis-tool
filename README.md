# Transmission Line Analysis & Modeling Tool

A complete Python-based tool for analyzing overhead transmission lines using industry-standard electrical engineering models.  
This project calculates RLC parameters, ABCD matrices, sending-end and receiving-end quantities, and supports full geometry-based analysis.

---

## 🚀 Features

### 🔷 Geometry Mode
 Computes **Geometrical Mean Distance (Dm)** from conductor coordinates  
 Computes **Self-GMD (Ds)** from conductor radius  
 Supports **bundle conductors**  
 Auto-calculates inductance (L) and capacitance (C)  
 Real tower geometry → real electrical parameters


🔷 Power Flow / Analysis Mode

Calculates series impedance (R + jX) and shunt admittance (jB) matrices

Computes ABCD parameters for short, medium, and long transmission lines

Calculates sending-end and receiving-end voltages and currents

Supports load flow calculations for balanced three-phase systems

Provides voltage regulation, line losses, and efficiency


🔷 User-Friendly Features

Interactive console-based interface

Accepts input in metric or imperial units

Provides clear, step-by-step output

Option to export results to CSV or Excel for documentation


🔷 Advanced Capabilities

Handles single-phase and three-phase lines

Accounts for conductor bundling and spacing variations

Incorporates Earth return effects in impedance calculations

Optional temperature correction for resistance values

📚 Technologies Used

Python 3.13

NumPy for matrix and numerical operations

Pandas for data handling and exporting

Matplotlib (optional) for visualization of line parameters
