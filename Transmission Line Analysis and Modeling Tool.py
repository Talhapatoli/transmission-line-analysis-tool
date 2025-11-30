import math
import cmath
import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
from datetime import datetime

# ===============================
#  AAC Conductor Data Dictionary
# ===============================
aac_data = {
    "Waxwing": {"Aluminum_area_cmil": 266800, "Stranding": "18/1", "Layers_Al": 2, "Outside_diameter_in": 0.609, "GMR_ft": 0.0093},
    "Partridge": {"Aluminum_area_cmil": 266800, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.642, "GMR_ft": 0.0095},
    "Ostrich": {"Aluminum_area_cmil": 300000, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.680, "GMR_ft": 0.0097},
    "Merlin": {"Aluminum_area_cmil": 336400, "Stranding": "18/1", "Layers_Al": 2, "Outside_diameter_in": 0.684, "GMR_ft": 0.0097},
    "Linnet": {"Aluminum_area_cmil": 336400, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.721, "GMR_ft": 0.0099},
    "Oriole": {"Aluminum_area_cmil": 336400, "Stranding": "30/7", "Layers_Al": 2, "Outside_diameter_in": 0.741, "GMR_ft": 0.0100},
    "Chickadee": {"Aluminum_area_cmil": 397500, "Stranding": "18/1", "Layers_Al": 2, "Outside_diameter_in": 0.743, "GMR_ft": 0.0101},
    "Ibis": {"Aluminum_area_cmil": 397500, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.783, "GMR_ft": 0.0102},
    "Pelican": {"Aluminum_area_cmil": 477000, "Stranding": "18/1", "Layers_Al": 2, "Outside_diameter_in": 0.814, "GMR_ft": 0.0104},
    "Flicker": {"Aluminum_area_cmil": 477000, "Stranding": "24/7", "Layers_Al": 2, "Outside_diameter_in": 0.846, "GMR_ft": 0.0106},
    "Hawk": {"Aluminum_area_cmil": 477000, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.858, "GMR_ft": 0.0107},
    "Hen": {"Aluminum_area_cmil": 477000, "Stranding": "30/7", "Layers_Al": 2, "Outside_diameter_in": 0.883, "GMR_ft": 0.0108},
    "Osprey": {"Aluminum_area_cmil": 556500, "Stranding": "18/1", "Layers_Al": 2, "Outside_diameter_in": 0.879, "GMR_ft": 0.0109},
    "Parakeet": {"Aluminum_area_cmil": 556500, "Stranding": "24/7", "Layers_Al": 2, "Outside_diameter_in": 0.914, "GMR_ft": 0.0110},
    "Dove": {"Aluminum_area_cmil": 556500, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.927, "GMR_ft": 0.0111},
    "Rook": {"Aluminum_area_cmil": 636000, "Stranding": "24/7", "Layers_Al": 2, "Outside_diameter_in": 0.977, "GMR_ft": 0.0113},
    "Grosbeak": {"Aluminum_area_cmil": 636000, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 0.990, "GMR_ft": 0.0114},
    "Drake": {"Aluminum_area_cmil": 795000, "Stranding": "26/7", "Layers_Al": 2, "Outside_diameter_in": 1.108, "GMR_ft": 0.0118},
    "Tern": {"Aluminum_area_cmil": 795000, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.063, "GMR_ft": 0.0116},
    "Rail": {"Aluminum_area_cmil": 954000, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.165, "GMR_ft": 0.0119},
    "Cardinal": {"Aluminum_area_cmil": 954000, "Stranding": "54/7", "Layers_Al": 3, "Outside_diameter_in": 1.196, "GMR_ft": 0.0120},
    "Ortolan": {"Aluminum_area_cmil": 1033500, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.213, "GMR_ft": 0.0122},
    "Bluejay": {"Aluminum_area_cmil": 1113000, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.259, "GMR_ft": 0.0123},
    "Finch": {"Aluminum_area_cmil": 1113000, "Stranding": "54/19", "Layers_Al": 3, "Outside_diameter_in": 1.293, "GMR_ft": 0.0124},
    "Bittern": {"Aluminum_area_cmil": 1272000, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.345, "GMR_ft": 0.0126},
    "Pheasant": {"Aluminum_area_cmil": 1272000, "Stranding": "54/19", "Layers_Al": 3, "Outside_diameter_in": 1.382, "GMR_ft": 0.0127},
    "Bobolink": {"Aluminum_area_cmil": 1431000, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.427, "GMR_ft": 0.0129},
    "Plover": {"Aluminum_area_cmil": 1431000, "Stranding": "54/19", "Layers_Al": 3, "Outside_diameter_in": 1.465, "GMR_ft": 0.0130},
    "Lapwing": {"Aluminum_area_cmil": 1590000, "Stranding": "45/7", "Layers_Al": 3, "Outside_diameter_in": 1.502, "GMR_ft": 0.0132},
    "Falcon": {"Aluminum_area_cmil": 1590000, "Stranding": "54/19", "Layers_Al": 3, "Outside_diameter_in": 1.545, "GMR_ft": 0.0133},
    "Bluebird": {"Aluminum_area_cmil": 2156000, "Stranding": "84/19", "Layers_Al": 4, "Outside_diameter_in": 1.762, "GMR_ft": 0.0138}
}


# ==============================
# === UNIT CONVERSION SYSTEM ===
# ==============================
length_units = {
    "meters": 1.0,
    "centimeters": 0.01,
    "millimeters": 0.001,
    "feet": 0.3048,
    "inches": 0.0254,
    "miles": 1609.344,
    "kilometers": 1000.0
}

# ==============================
# === CORE CALCULATION ENGINE ===
# ==============================

def feet_to_meters(feet):
    return feet * 0.3048

def inches_to_meters(inches):
    return inches * 0.0254

def convert_length(value, from_unit, to_unit):
    """Convert length between different units"""
    base_value = value * length_units[from_unit]
    return base_value / length_units[to_unit]

def calc_Ds(n, r_prime_ft, d_m):
    """Calculate equivalent GMR for bundled conductors"""
    correction_factors = {1: 1.0, 2: 1.0, 3: 1.0, 4: 1.09, 5: 1.178, 6: 1.25, 7: 1.32, 8: 1.38}
    k = correction_factors.get(n, 1.0)
    
    r_prime_m = feet_to_meters(r_prime_ft)
    
    if n == 1:
        return r_prime_m
    else:
        return k * ((r_prime_m * (d_m ** (n - 1))) ** (1 / n))

def calc_RLC_parameters(freq, Ds_m, Dm_m, conductor_area_cmil=None, temperature=20, output_units="per_km"):
    """Calculate R, XL, XC with multiple output unit options"""
    if conductor_area_cmil:
        area_m2 = conductor_area_cmil * 5.067e-10
        rho_al = 2.82e-8
        R_per_m = rho_al / area_m2
        alpha = 0.00403
        R_per_m *= (1 + alpha * (temperature - 20))
    else:
        R_per_m = 0.1e-3
    
    # Calculate base parameters per meter
    R_per_meter = R_per_m
    XL_per_meter = 2 * math.pi * freq * 2e-7 * math.log(Dm_m / Ds_m)
    C_per_m = (2 * math.pi * 8.854e-12) / math.log(Dm_m / Ds_m)
    XC_per_meter = 1 / (2 * math.pi * freq * C_per_m)
    
    # Convert to desired output units
    unit_conversion = {
        "per_m": 1.0,
        "per_km": 1000.0,
        "per_mile": 1609.344,
        "per_100km": 100000.0,
        "per_100mile": 160934.4
    }
    
    conversion = unit_conversion.get(output_units, 1000.0)
    
    R_output = R_per_meter * conversion
    XL_output = XL_per_meter * conversion
    XC_output = XC_per_meter / conversion * 1e-6
    
    return R_output, XL_output, XC_output

def calc_RLC_per_km(freq, Ds_m, Dm_m, conductor_area_cmil=None, temperature=20):
    return calc_RLC_parameters(freq, Ds_m, Dm_m, conductor_area_cmil, temperature, "per_km")

def determine_line_model(length_km, freq):
    """Intelligently determine the appropriate line model"""
    if length_km < 80:
        return "short", "Short Line Model (l < 80 km)"
    elif length_km <= 240:
        return "pi", f"Medium Line Model (80-240 km) - π-model at {freq} Hz"
    else:
        return "long", f"Long Line Model (l > 240 km) - Distributed parameters at {freq} Hz"

def calc_ABCD_parameters(length_km, R_per_km, XL_per_km, XC_MΩ_km, freq, model_type):
    """Calculate ABCD parameters with numerical safeguards"""
    ω = 2 * math.pi * freq
    Z_per_km = complex(R_per_km, XL_per_km)
    
    # Safeguard for capacitive reactance to prevent division by zero
    if XC_MΩ_km <= 0 or XC_MΩ_km > 1e6:
        Y_per_km = complex(0, 1e-12)  # Very small admittance for stability
    else:
        Y_per_km = complex(0, 1 / (XC_MΩ_km * 1e6))
    
    Z_total = Z_per_km * length_km
    Y_total = Y_per_km * length_km
    
    # Force short-line model for very short lines to avoid numerical issues
    if length_km < 10:
        model_type = "short"
    
    if model_type == "short":
        A = D = complex(1, 0)
        B = Z_total
        C = complex(0, 0)
        
    elif model_type == "pi":
        A = D = 1 + (Y_total * Z_total) / 2
        B = Z_total
        C = Y_total * (1 + (Y_total * Z_total) / 4)
        
    else:  # long line
        try:
            γ = cmath.sqrt(Z_total * Y_total)
            Zc = cmath.sqrt(Z_total / Y_total)
            l_m = length_km * 1000
            
            # Check for numerical stability before using hyperbolic functions
            if abs(γ * l_m) > 10:  # Too large for stable hyperbolic calculation
                # Fallback to medium line model for numerical stability
                A = D = 1 + (Y_total * Z_total) / 2
                B = Z_total
                C = Y_total * (1 + (Y_total * Z_total) / 4)
            else:
                A = D = cmath.cosh(γ * l_m)
                B = Zc * cmath.sinh(γ * l_m)
                C = (1 / Zc) * cmath.sinh(γ * l_m)
        except:
            # Fallback to medium line if long line calculation fails
            A = D = 1 + (Y_total * Z_total) / 2
            B = Z_total
            C = Y_total * (1 + (Y_total * Z_total) / 4)
    
    return A, B, C, D

def calc_power_flow(Vr_kV, Pr_MW, pf, A, B, C, D, lead_lag="lag"):
    """Calculate power flow with improved numerical stability"""
    Vr_phase = Vr_kV * 1000 / math.sqrt(3)
    
    # Handle very light loads better to avoid division issues
    if Pr_MW < 0.1:
        Ir_mag = 1.0  # Small current for numerical stability
    else:
        Ir_mag = (Pr_MW * 1e6) / (math.sqrt(3) * Vr_kV * 1000 * pf)
    
    θ = math.acos(pf)
    if lead_lag == "lead":
        θ = -θ
    
    Ir = cmath.rect(Ir_mag, -θ)
    Vs = A * Vr_phase + B * Ir
    Is = C * Vr_phase + D * Ir
    
    Vs_line_kV = abs(Vs) * math.sqrt(3) / 1000
    Is_line_A = abs(Is)
    
    Ss = 3 * Vs * Is.conjugate()
    Ps = Ss.real / 1e6
    Qs = Ss.imag / 1e6
    
    # Improved efficiency calculation for edge cases
    if Pr_MW < 0.1:
        efficiency = 99.9  # Assume high efficiency for very light loads
    else:
        efficiency = (Pr_MW / Ps * 100) if Ps > 0 and Pr_MW > 0 else 99.9
    
    # Safeguard for voltage regulation calculation
    Vr_no_load = abs(Vs) / abs(A) if abs(A) > 1e-12 else abs(Vs)
    voltage_regulation = ((Vr_no_load - abs(Vr_phase)) / abs(Vr_phase) * 100)
    
    return {
        'Vs_kV': Vs_line_kV,
        'Is_A': Is_line_A,
        'Ps_MW': Ps,
        'Qs_MVAR': Qs,
        'efficiency': efficiency,
        'voltage_regulation': voltage_regulation,
        'Vs_phase': Vs,
        'Is_phase': Is,
        'Vr_phase': Vr_phase,
        'Ir_phase': Ir
    }

def plot_unified_phasor_diagram(Vs_phase, Is_phase, Vr_phase, Ir_phase):
    """Plot all 4 phasors in one professional diagram"""
    if Vs_phase is None or Is_phase is None:
        print("No valid phasor data available")
        return
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Normalize for better visualization
    max_voltage = max(abs(Vs_phase), abs(Vr_phase)) if max(abs(Vs_phase), abs(Vr_phase)) > 0 else 1
    max_current = max(abs(Is_phase), abs(Ir_phase)) if max(abs(Is_phase), abs(Ir_phase)) > 0 else 1
    
    scale_voltage = 1.0 / max_voltage
    scale_current = 1.0 / max_current * 0.7
    
    # Scale phasors
    Vs_scaled = Vs_phase * scale_voltage
    Vr_scaled = Vr_phase * scale_voltage
    Is_scaled = Is_phase * scale_current
    Ir_scaled = Ir_phase * scale_current
    
    # Colors for different phasors
    colors = {'Vs': '#FF4444', 'Is': '#FFAA00', 'Vr': '#4444FF', 'Ir': '#00AA00'}
    
    # Calculate angles in degrees
    vs_angle = math.degrees(cmath.phase(Vs_phase))
    is_angle = math.degrees(cmath.phase(Is_phase))
    vr_angle = math.degrees(cmath.phase(Vr_phase))
    ir_angle = math.degrees(cmath.phase(Ir_phase))
    
    # Plot voltage phasors
    ax.quiver(0, 0, Vs_scaled.real, Vs_scaled.imag, 
              angles='xy', scale_units='xy', scale=1, 
              color=colors['Vs'], width=0.018, linewidth=2.5,
              label=f'Vs ({abs(Vs_phase)/1000:.1f} kV, ∠{vs_angle:.1f}°)')
    
    ax.quiver(0, 0, Vr_scaled.real, Vr_scaled.imag,
              angles='xy', scale_units='xy', scale=1,
              color=colors['Vr'], width=0.018, linewidth=2.5,
              label=f'Vr ({abs(Vr_phase)/1000:.1f} kV, ∠{vr_angle:.1f}°)')
    
    # Plot current phasors
    ax.quiver(0, 0, Is_scaled.real, Is_scaled.imag,
              angles='xy', scale_units='xy', scale=1,
              color=colors['Is'], width=0.012, linewidth=2.0,
              label=f'Is ({abs(Is_phase):.1f} A, ∠{is_angle:.1f}°)')
    
    ax.quiver(0, 0, Ir_scaled.real, Ir_scaled.imag,
              angles='xy', scale_units='xy', scale=1,
              color=colors['Ir'], width=0.012, linewidth=2.0,
              label=f'Ir ({abs(Ir_phase):.1f} A, ∠{ir_angle:.1f}°)')
    
    # Set plot limits
    max_magnitude = max(abs(Vs_scaled), abs(Vr_scaled), abs(Is_scaled), abs(Ir_scaled))
    limit = max_magnitude * 1.3
    
    # Create circular axes
    circle_radius = max_magnitude * 1.1
    circle = plt.Circle((0, 0), circle_radius, fill=False, color='gray', 
                       linestyle='--', alpha=0.4, linewidth=1)
    ax.add_patch(circle)
    
    # Add angle markers
    angles_deg = [0, 30, 45, 60, 90, 120, 135, 150, 180, 210, 225, 240, 270, 300, 315, 330]
    for angle in angles_deg:
        angle_rad = math.radians(angle)
        x_text = circle_radius * 1.08 * math.cos(angle_rad)
        y_text = circle_radius * 1.08 * math.sin(angle_rad)
        ax.text(x_text, y_text, f'{angle}°', ha='center', va='center', 
                fontsize=8, alpha=0.7, color='darkblue')
    
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    
    # Professional styling
    ax.grid(True, linestyle='--', alpha=0.4, linewidth=0.5)
    ax.set_aspect('equal')
    ax.axhline(y=0, color='k', linewidth=1, alpha=0.5)
    ax.axvline(x=0, color='k', linewidth=1, alpha=0.5)
    
    ax.set_title('UNIFIED PHASOR DIAGRAM - TRANSMISSION LINE ANALYSIS', 
                 fontsize=14, fontweight='bold', pad=20, color='darkblue')
    ax.set_xlabel('Real Component (Reference Axis)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Imaginary Component', fontsize=11, fontweight='bold')
    
    ax.legend(loc='upper right', fontsize=10, framealpha=0.95, 
              edgecolor='black', fancybox=True)
    
    plt.tight_layout()
    plt.show()

# ==============================
# === GEOMETRY DESIGNER CLASS ===
# ==============================

class TransmissionLineGeometryDesigner:
    def __init__(self, parent, main_app):
        self.parent = parent
        self.main_app = main_app
        self.cell_size = 8  # pixels per meter
        self.grid_size = 200  # 200x200 meter grid
        self.canvas_size = self.grid_size * self.cell_size
        
        # 1 grid unit = 1 meter exactly
        self.grid_to_meters_scale = 1.0
        
        # Core data structure - store decimal coordinates
        self.circuits = {}
        self.current_circuit = 1
        self.current_phase = 'A'
        self.symbol_mode = 'symbols'
        
        self.setup_default_circuits()
        self.setup_simplified_ui()
        self.draw_grid()
        
    def setup_default_circuits(self):
        """Create 3 default circuits"""
        for i in range(1, 4):
            self.circuits[i] = {
                'name': f'Circuit {i}',
                'A': [], 'B': [], 'C': [], 'G': [],
                'style': 'solid'
            }

    def setup_simplified_ui(self):
        """TLGD with visual controls AND coordinate input"""
        self.main_frame = ttk.Frame(self.parent)
        self.main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Left panel with grid
        left_panel = ttk.Frame(self.main_frame)
        left_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.setup_controls(left_panel)
        self.setup_grid_canvas(left_panel)
        
        # Right panel with coordinate input
        right_panel = ttk.Frame(self.main_frame, width=300)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, padx=(10, 0))
        right_panel.pack_propagate(False)
        
        self.setup_coordinate_input(right_panel)
        
    def setup_coordinate_input(self, parent):
        """Setup coordinate input for precise decimal placement"""
        coord_frame = ttk.LabelFrame(parent, text="Manual Placement (Meters)", padding=10)
        coord_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Instructions
        ttk.Label(coord_frame, text="Enter coordinates in meters (0-199.9):", 
                 font=("Arial", 9, "bold")).pack(anchor='w', pady=(0, 5))
        
        # Coordinate inputs with decimal support
        input_frame = ttk.Frame(coord_frame)
        input_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(input_frame, text="X (m):").grid(row=0, column=0, padx=2)
        self.x_var = tk.StringVar()
        self.x_entry = ttk.Entry(input_frame, textvariable=self.x_var, width=10)
        self.x_entry.grid(row=0, column=1, padx=2)
        
        ttk.Label(input_frame, text="Y (m):").grid(row=0, column=2, padx=2)
        self.y_var = tk.StringVar()
        self.y_entry = ttk.Entry(input_frame, textvariable=self.y_var, width=10)
        self.y_entry.grid(row=0, column=3, padx=2)
        
        # Action buttons
        button_frame = ttk.Frame(coord_frame)
        button_frame.pack(fill=tk.X, pady=5)
        
        ttk.Button(button_frame, text="Place Conductor", 
                  command=self.place_at_coordinates).pack(side=tk.LEFT, padx=2)
        ttk.Button(button_frame, text="Remove Conductor", 
                  command=self.remove_at_coordinates).pack(side=tk.LEFT, padx=2)
        
        # Bind Enter key to place conductor
        self.x_entry.bind('<Return>', lambda e: self.place_at_coordinates())
        self.y_entry.bind('<Return>', lambda e: self.place_at_coordinates())
        
        # Status for coordinate operations
        self.coord_status = tk.StringVar(value="Enter coordinates in meters and click Place")
        ttk.Label(coord_frame, textvariable=self.coord_status, 
                 font=("Arial", 8)).pack(anchor='w')
        
        # Quick placement guide with realistic transmission line spacing
        guide_frame = ttk.LabelFrame(coord_frame, text="Typical Transmission Line Spacing", padding=5)
        guide_frame.pack(fill=tk.X, pady=(10, 0))
        
        guide_text = """Common configurations (in meters):
• 3-phase horizontal: (5,95), (10,95), (15,95)  [5m spacing]
• 3-phase vertical: (10,5), (10,10), (10,15)    [5m spacing]  
• 3-phase delta: (8,8), (12,8), (10,12)         [~4-5m spacing]
• Bundle conductors: Multiple at same phase position"""
        
        guide_label = ttk.Label(guide_frame, text=guide_text, 
                               font=("Arial", 8), justify=tk.LEFT)
        guide_label.pack(anchor='w')

    def place_at_coordinates(self):
        """Place conductor at specified decimal coordinates"""
        try:
            x = float(self.x_var.get())
            y = float(self.y_var.get())
            if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
                success = self.place_conductor(x, y)
                if success:
                    self.draw_conductors()
                    self.update_status(f"Placed {self.get_current_phase_name()} at ({x:.2f}m, {y:.2f}m)")
                    self.coord_status.set(f"Placed at ({x:.2f}m, {y:.2f}m)")
                    
                    # Clear inputs after successful placement
                    self.x_var.set("")
                    self.y_var.set("")
                    self.x_entry.focus()
                else:
                    self.coord_status.set("Placement failed - position occupied?")
            else:
                self.coord_status.set(f"Coordinates must be 0-{self.grid_size-1} meters")
        except ValueError:
            self.coord_status.set("Enter valid numbers (decimals allowed)")

    def remove_at_coordinates(self):
        """Remove conductor at specified decimal coordinates"""
        try:
            x = float(self.x_var.get())
            y = float(self.y_var.get())
            if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
                success = self.remove_conductor(x, y)
                if success:
                    self.draw_conductors()
                    self.update_status(f"Removed from ({x:.2f}m, {y:.2f}m)")
                    self.coord_status.set(f"Removed from ({x:.2f}m, {y:.2f}m)")
                    
                    # Clear inputs
                    self.x_var.set("")
                    self.y_var.set("")
                else:
                    self.coord_status.set("No conductor found at specified coordinates")
            else:
                self.coord_status.set(f"Coordinates must be 0-{self.grid_size-1} meters")
        except ValueError:
            self.coord_status.set("Enter valid numbers (decimals allowed)")

    def setup_controls(self, parent):
        """Only visual controls - no electrical parameters"""
        self.control_frame = ttk.Frame(parent)
        self.control_frame.pack(fill=tk.X, pady=(0, 10))
        
        # Phase selection
        ttk.Label(self.control_frame, text="Current Phase:").grid(row=0, column=0, padx=5)
        self.phase_var = tk.StringVar(value='A (Red)')
        phase_combo = ttk.Combobox(self.control_frame, textvariable=self.phase_var, 
                                  values=['A (Red)', 'B (Yellow)', 'C (Blue)', 'G (Ground)'],
                                  width=12, state='readonly')
        phase_combo.grid(row=0, column=1, padx=5)
        phase_combo.bind('<<ComboboxSelected>>', self.on_phase_change)
        
        # Circuit selection
        ttk.Label(self.control_frame, text="Circuit:").grid(row=0, column=2, padx=5)
        self.circuit_var = tk.StringVar()
        self.circuit_combo = ttk.Combobox(self.control_frame, textvariable=self.circuit_var,
                                         width=12, state='readonly')
        self.circuit_combo.grid(row=0, column=3, padx=5)
        self.circuit_combo.bind('<<ComboboxSelected>>', self.on_circuit_change)
        
        # Circuit management
        ttk.Button(self.control_frame, text="+ Circuit", 
                  command=self.add_circuit).grid(row=0, column=4, padx=2)
        ttk.Button(self.control_frame, text="- Circuit", 
                  command=self.remove_circuit).grid(row=0, column=5, padx=2)
        
        # Display mode
        ttk.Label(self.control_frame, text="Display:").grid(row=1, column=0, padx=5, pady=5)
        self.symbol_var = tk.StringVar(value='symbols')
        ttk.Radiobutton(self.control_frame, text="Symbols", variable=self.symbol_var,
                       value='symbols', command=self.on_symbol_change).grid(row=1, column=1, padx=2)
        ttk.Radiobutton(self.control_frame, text="Letters", variable=self.symbol_var,
                       value='letters', command=self.on_symbol_change).grid(row=1, column=2, padx=2)
        
        # Action buttons
        ttk.Button(self.control_frame, text="Clear Phase", 
                  command=self.clear_current_phase).grid(row=1, column=3, padx=2)
        ttk.Button(self.control_frame, text="Clear All", 
                  command=self.clear_all_phases).grid(row=1, column=4, padx=2)
        ttk.Button(self.control_frame, text="Reset Grid", 
                  command=self.reset_grid).grid(row=1, column=5, padx=2)
        
        # Transfer button
        ttk.Button(self.control_frame, text="Calculate Dm & Transfer", 
                  command=self.calculate_and_transfer).grid(row=0, column=6, padx=10)
        
        # Status
        self.status_var = tk.StringVar(value="Place conductors then click 'Calculate Dm & Transfer'")
        ttk.Label(self.control_frame, textvariable=self.status_var).grid(row=2, column=0, columnspan=7, pady=5, sticky='w')
        
        self.update_circuit_list()

    def setup_grid_canvas(self, parent):
        """Setup grid canvas"""
        self.canvas_frame = ttk.Frame(parent)
        self.canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.h_scroll = ttk.Scrollbar(self.canvas_frame, orient=tk.HORIZONTAL)
        self.h_scroll.pack(side=tk.BOTTOM, fill=tk.X)
        self.v_scroll = ttk.Scrollbar(self.canvas_frame, orient=tk.VERTICAL)
        self.v_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.canvas = tk.Canvas(
            self.canvas_frame,
            width=800,
            height=600,
            bg='white',
            xscrollcommand=self.h_scroll.set,
            yscrollcommand=self.v_scroll.set,
            scrollregion=(0, 0, self.canvas_size, self.canvas_size)
        )
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.h_scroll.config(command=self.canvas.xview)
        self.v_scroll.config(command=self.canvas.yview)
        
        self.canvas.bind("<Button-1>", self.on_canvas_click)
        self.canvas.bind("<Button-3>", self.on_canvas_right_click)
        self.canvas.bind("<Motion>", self.on_canvas_motion)

    def calculate_and_transfer(self):
        """TLGD only calculates Dm and transfers to main app"""
        try:
            # Calculate ONLY Dm from geometry
            dm = self.calculate_geometric_mean_distance()
            
            # Transfer to main app
            self.main_app.receive_geometry_parameters({
                'dm': dm,
                'source': 'geometry_designer'
            })
            
            self.update_status(f"Dm = {dm:.3f}m transferred to main calculator")
            
        except Exception as e:
            self.update_status(f"Error calculating Dm: {str(e)}")

    def calculate_geometric_mean_distance(self) -> float:
        circuit = self.circuits[self.current_circuit]
        centroids = {}
        for phase in ['A', 'B', 'C']:
            if circuit[phase]:
                # Use decimal coordinates directly (already in meters)
                x_avg = sum(pos[0] for pos in circuit[phase]) / len(circuit[phase])
                y_avg = sum(pos[1] for pos in circuit[phase]) / len(circuit[phase])
                centroids[phase] = (x_avg, y_avg)
        
        distances = []
        phases = list(centroids.keys())
        for i in range(len(phases)):
            for j in range(i + 1, len(phases)):
                x1, y1 = centroids[phases[i]]
                x2, y2 = centroids[phases[j]]
                distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                distances.append(distance)
        
        if distances and all(d > 0 for d in distances):
            product = 1.0
            for dist in distances:
                product *= dist
            dm = product ** (1 / len(distances))
            return max(dm, 0.1)
        else:
            return 5.0

    def draw_grid(self):
        self.canvas.delete("grid")
        for i in range(0, self.grid_size + 1):
            x = i * self.cell_size
            self.canvas.create_line(x, 0, x, self.canvas_size, fill='#E0E0E0', tags="grid", width=1)
            y = i * self.cell_size
            self.canvas.create_line(0, y, self.canvas_size, y, fill='#E0E0E0', tags="grid", width=1)
        
        # Label every 10 meters
        for i in range(0, self.grid_size + 1, 10):
            x = i * self.cell_size
            self.canvas.create_line(x, 0, x, self.canvas_size, fill='#B0B0B0', tags="grid", width=2)
            y = i * self.cell_size  
            self.canvas.create_line(0, y, self.canvas_size, y, fill='#B0B0B0', tags="grid", width=2)
            
            if i > 0 and i < self.grid_size:
                self.canvas.create_text(x, 10, text=f"{i}m", fill='#666666', tags="grid", font=('Arial', 8))
                self.canvas.create_text(10, y, text=f"{i}m", fill='#666666', tags="grid", font=('Arial', 8))
        
        self.draw_conductors()

    def draw_conductors(self):
        self.canvas.delete("conductor")
        circuit = self.circuits[self.current_circuit]
        colors = {'A': '#FF4444', 'B': '#FFAA00', 'C': '#4444FF', 'G': '#666666'}
        
        for phase in ['A', 'B', 'C', 'G']:
            for x, y in circuit[phase]:
                # Convert meter coordinates to canvas pixels (1:1 scale)
                canvas_x = x * self.cell_size
                canvas_y = y * self.cell_size
                radius = self.cell_size // 2 - 1
                color = colors[phase]
                
                if self.symbol_mode == 'symbols':
                    if phase == 'A':
                        self.canvas.create_oval(canvas_x - radius, canvas_y - radius, canvas_x + radius, canvas_y + radius, fill=color, outline='black', width=1, tags="conductor")
                    elif phase == 'B':
                        points = [canvas_x, canvas_y - radius, canvas_x - radius, canvas_y + radius, canvas_x + radius, canvas_y + radius]
                        self.canvas.create_polygon(points, fill=color, outline='black', width=1, tags="conductor")
                    elif phase == 'C':
                        self.canvas.create_rectangle(canvas_x - radius, canvas_y - radius, canvas_x + radius, canvas_y + radius, fill=color, outline='black', width=1, tags="conductor")
                    else:
                        points = [canvas_x, canvas_y - radius, canvas_x + radius, canvas_y, canvas_x, canvas_y + radius, canvas_x - radius, canvas_y]
                        self.canvas.create_polygon(points, fill=color, outline='black', width=1, tags="conductor")
                else:
                    self.canvas.create_rectangle(canvas_x - radius, canvas_y - radius, canvas_x + radius, canvas_y + radius, fill=color, outline='black', width=1, tags="conductor")
                    self.canvas.create_text(canvas_x, canvas_y, text=phase, fill='white', font=('Arial', 8, 'bold'), tags="conductor")

    def place_conductor(self, x: float, y: float) -> bool:
        """Place conductor at decimal coordinates (meters)"""
        if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
            # Store as decimal coordinates
            self.circuits[self.current_circuit][self.current_phase].append((x, y))
            return True
        return False

    def remove_conductor(self, x: float, y: float) -> bool:
        """Remove conductor at decimal coordinates with tolerance"""
        conductors = self.circuits[self.current_circuit][self.current_phase]
        tolerance = 0.1  # 10cm tolerance for removal
        
        for i, (cx, cy) in enumerate(conductors):
            if abs(cx - x) < tolerance and abs(cy - y) < tolerance:
                conductors.pop(i)
                return True
        return False

    def on_canvas_click(self, event):
        """Convert click to decimal meter coordinates"""
        canvas_x = self.canvas.canvasx(event.x)
        canvas_y = self.canvas.canvasy(event.y)
        
        # Convert pixels to meters (1:1 scale)
        x = canvas_x / self.cell_size
        y = canvas_y / self.cell_size
        
        if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
            success = self.place_conductor(x, y)
            if success:
                self.draw_conductors()
                self.update_status(f"Placed {self.get_current_phase_name()} at ({x:.2f}m, {y:.2f}m)")

    def on_canvas_right_click(self, event):
        """Convert right-click to decimal meter coordinates for removal"""
        canvas_x = self.canvas.canvasx(event.x)
        canvas_y = self.canvas.canvasy(event.y)
        
        # Convert pixels to meters (1:1 scale)
        x = canvas_x / self.cell_size
        y = canvas_y / self.cell_size
        
        if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
            success = self.remove_conductor(x, y)
            if success:
                self.draw_conductors()
                self.update_status(f"Removed from ({x:.2f}m, {y:.2f}m)")

    def on_canvas_motion(self, event):
        """Show decimal meter coordinates during mouse movement"""
        canvas_x = self.canvas.canvasx(event.x)
        canvas_y = self.canvas.canvasy(event.y)
        
        x = canvas_x / self.cell_size
        y = canvas_y / self.cell_size
        
        if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
            self.update_status(f"Position: ({x:.2f}m, {y:.2f}m) - Click to place, Right-click to remove")

    def on_phase_change(self, event):
        phase_map = {'A (Red)': 'A', 'B (Yellow)': 'B', 'C (Blue)': 'C', 'G (Ground)': 'G'}
        selected = self.phase_var.get()
        if selected in phase_map:
            self.current_phase = phase_map[selected]
            self.update_status(f"Switched to {selected}")

    def on_circuit_change(self, event):
        """Handle circuit selection change"""
        selected = self.circuit_var.get()
        if selected and selected.startswith('Circuit '):
            try:
                circuit_id = int(selected.split(' ')[1])
                self.current_circuit = circuit_id
                self.draw_conductors()
                self.update_status(f"Switched to {selected}")
            except (ValueError, IndexError):
                pass

    def on_symbol_change(self):
        self.symbol_mode = self.symbol_var.get()
        self.draw_conductors()
        self.update_status(f"Switched to {self.symbol_mode} mode")

    def add_circuit(self):
        new_id = max(self.circuits.keys()) + 1 if self.circuits else 1
        self.circuits[new_id] = {
            'name': f'Circuit {new_id}',
            'A': [], 'B': [], 'C': [], 'G': [],
            'style': 'solid'
        }
        self.update_circuit_list()
        self.current_circuit = new_id
        self.circuit_var.set(f"Circuit {new_id}")
        self.draw_conductors()
        self.update_status(f"Added Circuit {new_id}")

    def remove_circuit(self):
        if len(self.circuits) > 1:
            del self.circuits[self.current_circuit]
            self.current_circuit = next(iter(self.circuits.keys()))
            self.update_circuit_list()
            self.circuit_var.set(f"Circuit {self.current_circuit}")
            self.draw_conductors()
            self.update_status("Removed circuit")

    def clear_current_phase(self):
        self.circuits[self.current_circuit][self.current_phase].clear()
        self.draw_conductors()
        self.update_status(f"Cleared {self.get_current_phase_name()}")

    def clear_all_phases(self):
        circuit = self.circuits[self.current_circuit]
        for phase in ['A', 'B', 'C', 'G']:
            circuit[phase].clear()
        self.draw_conductors()
        self.update_status("Cleared all phases")

    def reset_grid(self):
        for circuit in self.circuits.values():
            for phase in ['A', 'B', 'C', 'G']:
                circuit[phase].clear()
        self.draw_conductors()
        self.update_status("Grid reset")

    def update_circuit_list(self):
        circuits = [f"Circuit {cid}" for cid in sorted(self.circuits.keys())]
        self.circuit_combo['values'] = circuits
        if circuits and not self.circuit_var.get():
            self.circuit_var.set(f"Circuit {self.current_circuit}")

    def get_current_phase_name(self):
        names = {'A': 'Phase A (Red)', 'B': 'Phase B (Yellow)', 'C': 'Phase C (Blue)', 'G': 'Ground Wire'}
        return names.get(self.current_phase, 'Unknown Phase')

    def update_status(self, message: str):
        self.status_var.set(message)

    def update_display(self):
        self.update_circuit_list()
        self.draw_conductors()

# ==============================
# === MAIN APPLICATION CLASS ===
# ==============================

class TransmissionLineCalculator:
    def __init__(self, root):
        self.root = root
        self.root.title("Professional Transmission Line Calculator")
        self.root.geometry("1400x900")
        self.root.config(bg="#f0f0f0")
        
        self.last_calculation = None
        self.setup_gui()
        
    def setup_gui(self):
        # Main container
        main_container = tk.PanedWindow(self.root, orient=tk.HORIZONTAL, bg="#f0f0f0")
        main_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Left panel - Input controls
        left_frame = tk.Frame(main_container, bg="#f0f0f0", relief=tk.RAISED, bd=1)
        main_container.add(left_frame, stretch="always")
        
        # Header
        header = tk.Label(left_frame, text="Transmission Line Analysis System", 
                         font=("Arial", 16, "bold"), bg="#2c3e50", fg="white", pady=15)
        header.pack(fill=tk.X, pady=(0, 10))
        
        # Mode selection
        mode_frame = tk.Frame(left_frame, bg="#f0f0f0")
        mode_frame.pack(fill=tk.X, pady=10, padx=10)
        
        tk.Label(mode_frame, text="Analysis Mode:", font=("Arial", 12, "bold"), 
                bg="#f0f0f0").pack(side=tk.LEFT, padx=5)
        
        self.mode_var = tk.StringVar(value="geometry")
        tk.Radiobutton(mode_frame, text="Geometry Mode", variable=self.mode_var, 
                      value="geometry", command=self.switch_mode, bg="#f0f0f0",
                      font=("Arial", 10)).pack(side=tk.LEFT, padx=5)
        tk.Radiobutton(mode_frame, text="Power Flow Mode", variable=self.mode_var, 
                      value="powerflow", command=self.switch_mode, bg="#f0f0f0",
                      font=("Arial", 10)).pack(side=tk.LEFT, padx=5)
        
        # Input frame with scrollbar
        input_container = tk.Frame(left_frame, bg="#f0f0f0")
        input_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        input_canvas = tk.Canvas(input_container, bg="#f0f0f0", highlightthickness=0)
        scrollbar = tk.Scrollbar(input_container, orient="vertical", command=input_canvas.yview)
        scrollable_frame = tk.Frame(input_canvas, bg="#f0f0f0")
        
        scrollable_frame.bind("<Configure>", lambda e: input_canvas.configure(scrollregion=input_canvas.bbox("all")))
        input_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        input_canvas.configure(yscrollcommand=scrollbar.set)
        
        input_canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        self.input_frame = scrollable_frame
        
        # Control buttons
        button_frame = tk.Frame(left_frame, bg="#f0f0f0")
        button_frame.pack(fill=tk.X, pady=10, padx=10)
        
        self.calc_btn = tk.Button(button_frame, text="Calculate", command=self.calculate,
                                 font=("Arial", 12, "bold"), bg="#27ae60", fg="white", width=12)
        self.calc_btn.pack(side=tk.LEFT, padx=5)
        
        self.clear_btn = tk.Button(button_frame, text="Clear", command=self.clear_results,
                                 font=("Arial", 12), bg="#e74c3c", fg="white", width=10)
        self.clear_btn.pack(side=tk.LEFT, padx=5)
        
        self.save_btn = tk.Button(button_frame, text="Save Results", command=self.save_results,
                                 font=("Arial", 12), bg="#3498db", fg="white", width=12)
        self.save_btn.pack(side=tk.LEFT, padx=5)
        
        self.phasor_btn = tk.Button(button_frame, text="Unified Phasor", command=self.plot_unified_phasors,
                                   font=("Arial", 12), bg="#9b59b6", fg="white", width=12)
        self.phasor_btn.pack(side=tk.LEFT, padx=5)
        
        self.geometry_btn = tk.Button(button_frame, text="Plot Geometry", 
                                     command=self.plot_transmission_geometry,
                                     font=("Arial", 12), bg="#e67e22", fg="white", width=12)
        self.geometry_btn.pack(side=tk.LEFT, padx=5)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready - Select mode and enter parameters")
        status_bar = tk.Label(left_frame, textvariable=self.status_var, 
                             font=("Arial", 9), bg="#bdc3c7", fg="#2c3e50", 
                             relief="sunken", bd=1)
        status_bar.pack(fill=tk.X, side=tk.BOTTOM, pady=(5, 0))
        
        # Right panel - Results & Plots
        right_panel = tk.PanedWindow(main_container, orient=tk.VERTICAL, bg="#f0f0f0")
        main_container.add(right_panel, stretch="always")
        
        # Results display
        results_frame = tk.Frame(right_panel, bg="#f0f0f0", relief=tk.RAISED, bd=1)
        right_panel.add(results_frame, stretch="always")
        
        results_label = tk.Label(results_frame, text="Analysis Results", 
                                font=("Arial", 12, "bold"), bg="#f0f0f0")
        results_label.pack(pady=(5, 0))
        
        # Results text area
        results_text_frame = tk.Frame(results_frame, bg="#f0f0f0")
        results_text_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.result_text = tk.Text(results_text_frame, height=15, width=80, 
                                 font=("Consolas", 9), bg="white", relief="solid", bd=1)
        
        result_scrollbar = tk.Scrollbar(results_text_frame, command=self.result_text.yview)
        self.result_text.configure(yscrollcommand=result_scrollbar.set)
        
        self.result_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        result_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.build_geometry_inputs()
        
        # Set initial weights
        main_container.paneconfig(left_frame, minsize=400)
        right_panel.paneconfig(results_frame, minsize=300)

    def build_geometry_inputs(self):
        self.clear_input_frame()
        
        # Conductor selection
        tk.Label(self.input_frame, text="Conductor Type:", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=0, column=0, padx=5, pady=5, sticky="e")
        
        self.conductor_combo = ttk.Combobox(self.input_frame, 
                                          values=list(aac_data.keys()), 
                                          width=20, state="readonly")
        self.conductor_combo.grid(row=0, column=1, padx=5, pady=5)
        self.conductor_combo.bind('<<ComboboxSelected>>', self.on_conductor_select)
        
        # Bundle configuration
        tk.Label(self.input_frame, text="Subconductors per Bundle:", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=1, column=0, padx=5, pady=5, sticky="e")
        
        self.entry_n = tk.Entry(self.input_frame, width=10)
        self.entry_n.grid(row=1, column=1, padx=5, pady=5)
        self.entry_n.insert(0, "1")
        
        tk.Label(self.input_frame, text="Bundle Spacing:", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=2, column=0, padx=5, pady=5, sticky="e")
        
        self.entry_spacing = tk.Entry(self.input_frame, width=10)
        self.entry_spacing.grid(row=2, column=1, padx=5, pady=5)
        self.entry_spacing.insert(0, "0.3")
        
        self.spacing_unit = ttk.Combobox(self.input_frame, values=["meters"], 
                                       width=8, state="readonly")
        self.spacing_unit.grid(row=2, column=2, padx=5, pady=5)
        self.spacing_unit.set("meters")
        
        # Phase Spacing Dm
        dm_frame = tk.Frame(self.input_frame, bg="#f0f0f0")
        dm_frame.grid(row=3, column=0, columnspan=3, padx=5, pady=5, sticky="w")
        
        tk.Label(dm_frame, text="Phase Spacing (Dm):", bg="#f0f0f0",
                font=("Arial", 10)).pack(side=tk.LEFT, padx=5)
        
        self.entry_Dm = tk.Entry(dm_frame, width=10)
        self.entry_Dm.pack(side=tk.LEFT, padx=5)
        self.entry_Dm.insert(0, "5.0")
        
        # System parameters
        tk.Label(self.input_frame, text="Frequency (Hz):", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=4, column=0, padx=5, pady=5, sticky="e")
        
        self.entry_freq = tk.Entry(self.input_frame, width=10)
        self.entry_freq.grid(row=4, column=1, padx=5, pady=5)
        self.entry_freq.insert(0, "60")
        
        tk.Label(self.input_frame, text="Line Length (km):", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=5, column=0, padx=5, pady=5, sticky="e")
        
        self.entry_length = tk.Entry(self.input_frame, width=10)
        self.entry_length.grid(row=5, column=1, padx=5, pady=5)
        self.entry_length.insert(0, "100")

        # Parameter units
        tk.Label(self.input_frame, text="Parameter Units:", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=6, column=0, padx=5, pady=5, sticky="e")
        
        self.param_units = ttk.Combobox(self.input_frame, 
                                      values=["per_km", "per_mile", "per_m", "per_100km", "per_100mile"], 
                                      width=10, state="readonly")
        self.param_units.grid(row=6, column=1, padx=5, pady=5)
        self.param_units.set("per_km")
        
        # Output display units
        tk.Label(self.input_frame, text="Display Units:", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=7, column=0, padx=5, pady=5, sticky="e")
        
        self.output_unit = ttk.Combobox(self.input_frame, values=list(length_units.keys()), 
                                      width=10, state="readonly")
        self.output_unit.grid(row=7, column=1, padx=5, pady=5)
        self.output_unit.set("meters")
        
        # Conductor info display
        self.conductor_info = tk.Label(self.input_frame, text="Select a conductor to see details",
                                      bg="#f0f0f0", font=("Arial", 9), wraplength=300)
        self.conductor_info.grid(row=0, column=3, rowspan=4, padx=10, pady=5, sticky="w")

    def build_powerflow_inputs(self):
        self.clear_input_frame()
    
        # End selection
        tk.Label(self.input_frame, text="Known End:", bg="#f0f0f0",
                font=("Arial", 10, "bold")).grid(row=0, column=0, padx=5, pady=5, sticky="e")
    
        self.end_var = tk.StringVar(value="receiving")
        tk.Radiobutton(self.input_frame, text="Receiving End", variable=self.end_var, 
                      value="receiving", bg="#f0f0f0").grid(row=0, column=1, sticky="w")
        tk.Radiobutton(self.input_frame, text="Sending End", variable=self.end_var, 
                      value="sending", bg="#f0f0f0").grid(row=0, column=2, sticky="w")
    
        # Line parameters
        tk.Label(self.input_frame, text="Line Parameters:", bg="#f0f0f0",
                font=("Arial", 10, "bold")).grid(row=1, column=0, padx=5, pady=5, sticky="e")
    
        self.param_type = tk.StringVar(value="RLC_per_unit")
        param_combo = ttk.Combobox(self.input_frame, textvariable=self.param_type,
                            values=["RLC_per_unit", "RLC_total", "YZ", "geometry"], width=12)
        param_combo.grid(row=1, column=1, padx=5, pady=5)
        param_combo.bind('<<ComboboxSelected>>', self.on_param_type_change)
    
        # Parameter units
        tk.Label(self.input_frame, text="Parameter Units:", bg="#f0f0f0",
                font=("Arial", 10)).grid(row=1, column=2, padx=5, pady=5, sticky="e")
    
        self.powerflow_units = ttk.Combobox(self.input_frame,
                                          values=["per_km", "per_mile", "per_m"], 
                                          width=10, state="readonly")
        self.powerflow_units.grid(row=1, column=3, padx=5, pady=5)
        self.powerflow_units.set("per_km")
    
        # Parameter inputs frame
        self.param_frame = tk.Frame(self.input_frame, bg="#f0f0f0")
        self.param_frame.grid(row=2, column=0, columnspan=4, padx=5, pady=5, sticky="w")
    
        # End parameters
        tk.Label(self.input_frame, text="End Parameters:", bg="#f0f0f0",
                font=("Arial", 10, "bold")).grid(row=3, column=0, padx=5, pady=5, sticky="e")
    
        self.build_rlc_params()
        self.build_end_params()

    def build_rlc_params(self):
        for widget in self.param_frame.winfo_children():
            widget.destroy()
        
        selected_units = self.powerflow_units.get()
    
        unit_labels = {
            "per_km": ("Ω/km", "Ω/km", "MΩ·km"),
            "per_mile": ("Ω/mile", "Ω/mile", "MΩ·mile"), 
            "per_m": ("Ω/m", "Ω/m", "MΩ·m")
        }
    
        r_unit, xl_unit, xc_unit = unit_labels.get(selected_units, ("Ω/km", "Ω/km", "MΩ·km"))

        if self.param_type.get() == "RLC_per_unit":
            params = [
                (f"Resistance ({r_unit}):", "0.05"), 
                (f"Inductive Reactance ({xl_unit}):", "0.5"), 
                (f"Capacitive Reactance ({xc_unit}):", "0.2")
            ]
        elif self.param_type.get() == "RLC_total":
            params = [
                ("Total Resistance R (Ω):", "5.0"), 
                ("Total Inductive Reactance XL (Ω):", "50.0"),
                ("Total Capacitive Reactance XC (MΩ):", "20.0")
            ]
        elif self.param_type.get() == "YZ":
            params = [
                ("Total Impedance Z (Ω):", "5+50j"), 
                ("Total Admittance Y (S):", "0+0.02j")
            ]
        else:  # geometry
            self.build_geometry_param_inputs()
            return
            
        for i, (label, default) in enumerate(params):
            tk.Label(self.param_frame, text=label, bg="#f0f0f0",
                    font=("Arial", 9)).grid(row=i, column=0, padx=5, pady=2, sticky="e")
            entry = tk.Entry(self.param_frame, width=12)
            entry.grid(row=i, column=1, padx=5, pady=2)
            entry.insert(0, default)
            setattr(self, f"param_entry_{i}", entry)
        
        if self.param_type.get() in ["RLC_per_unit", "RLC_total"]:
            tk.Label(self.param_frame, text="Line Length (km):", bg="#f0f0f0",
                    font=("Arial", 9)).grid(row=len(params), column=0, padx=5, pady=2, sticky="e")
            self.entry_line_length = tk.Entry(self.param_frame, width=12)
            self.entry_line_length.grid(row=len(params), column=1, padx=5, pady=2)
            self.entry_line_length.insert(0, "100")

    def build_geometry_param_inputs(self):
        for widget in self.param_frame.winfo_children():
            widget.destroy()
            
        tk.Label(self.param_frame, text="Conductor Type:", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=0, column=0, padx=5, pady=2, sticky="e")
        
        self.geom_conductor_combo = ttk.Combobox(self.param_frame, 
                                               values=list(aac_data.keys()), 
                                               width=15, state="readonly")
        self.geom_conductor_combo.grid(row=0, column=1, padx=5, pady=2)
        
        # Bundle configuration
        tk.Label(self.param_frame, text="Bundle Conductors:", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=1, column=0, padx=5, pady=2, sticky="e")
        self.geom_bundle = tk.Entry(self.param_frame, width=12)
        self.geom_bundle.grid(row=1, column=1, padx=5, pady=2)
        self.geom_bundle.insert(0, "1")
        
        tk.Label(self.param_frame, text="Bundle Spacing (m):", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=2, column=0, padx=5, pady=2, sticky="e")
        self.geom_spacing = tk.Entry(self.param_frame, width=12)
        self.geom_spacing.grid(row=2, column=1, padx=5, pady=2)
        self.geom_spacing.insert(0, "0.3")
        
        # Phase Spacing
        dm_frame = tk.Frame(self.param_frame, bg="#f0f0f0")
        dm_frame.grid(row=3, column=0, columnspan=2, padx=5, pady=2, sticky="w")
        
        tk.Label(dm_frame, text="Phase Spacing (Dm):", bg="#f0f0f0",
                font=("Arial", 9)).pack(side=tk.LEFT, padx=5)
        
        self.geom_Dm = tk.Entry(dm_frame, width=12)
        self.geom_Dm.pack(side=tk.LEFT, padx=5)
        self.geom_Dm.insert(0, "5.0")
        
        # Other geometry parameters
        tk.Label(self.param_frame, text="Line Length (km):", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=5, column=0, padx=5, pady=2, sticky="e")
        self.geom_length = tk.Entry(self.param_frame, width=12)
        self.geom_length.grid(row=5, column=1, padx=5, pady=2)
        self.geom_length.insert(0, "100")
        
        tk.Label(self.param_frame, text="Frequency (Hz):", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=6, column=0, padx=5, pady=2, sticky="e")
        self.geom_freq = tk.Entry(self.param_frame, width=12)
        self.geom_freq.grid(row=6, column=1, padx=5, pady=2)
        self.geom_freq.insert(0, "60")
    
    def build_end_params(self):
        end_frame = tk.Frame(self.input_frame, bg="#f0f0f0")
        end_frame.grid(row=4, column=0, columnspan=4, padx=5, pady=5, sticky="w")
        
        # Voltage
        tk.Label(end_frame, text="Voltage (kV):", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=0, column=0, padx=5, pady=2, sticky="e")
        self.entry_voltage = tk.Entry(end_frame, width=12)
        self.entry_voltage.grid(row=0, column=1, padx=5, pady=2)
        self.entry_voltage.insert(0, "230")
        
        # Power
        tk.Label(end_frame, text="Power (MW):", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=0, column=2, padx=5, pady=2, sticky="e")
        self.entry_power = tk.Entry(end_frame, width=12)
        self.entry_power.grid(row=0, column=3, padx=5, pady=2)
        self.entry_power.insert(0, "100")
        
        # Power Factor
        tk.Label(end_frame, text="Power Factor:", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=1, column=0, padx=5, pady=2, sticky="e")
        self.entry_pf = tk.Entry(end_frame, width=12)
        self.entry_pf.grid(row=1, column=1, padx=5, pady=2)
        self.entry_pf.insert(0, "0.95")
        
        # PF Type
        self.pf_var = tk.StringVar(value="lag")
        tk.Radiobutton(end_frame, text="Lagging", variable=self.pf_var, 
                      value="lag", bg="#f0f0f0").grid(row=1, column=2, sticky="w")
        tk.Radiobutton(end_frame, text="Leading", variable=self.pf_var, 
                      value="lead", bg="#f0f0f0").grid(row=1, column=3, sticky="w")
        
        # Frequency
        tk.Label(end_frame, text="Frequency (Hz):", bg="#f0f0f0",
                font=("Arial", 9)).grid(row=2, column=0, padx=5, pady=2, sticky="e")
        self.entry_freq_pf = tk.Entry(end_frame, width=12)
        self.entry_freq_pf.grid(row=2, column=1, padx=5, pady=2)
        self.entry_freq_pf.insert(0, "60")

    def on_param_type_change(self, event):
        self.build_rlc_params()

    def clear_input_frame(self):
        for widget in self.input_frame.winfo_children():
            widget.destroy()

    def switch_mode(self):
        if self.mode_var.get() == "geometry":
            self.build_geometry_inputs()
            self.status_var.set("Geometry Mode: Calculate parameters from conductor geometry")
        else:
            self.build_powerflow_inputs()
            self.status_var.set("Power Flow Mode: Analyze power flow with given parameters")
    
        self.phasor_btn.config(state="disabled")

    def on_conductor_select(self, event):
        conductor = self.conductor_combo.get()
        if conductor in aac_data:
            data = aac_data[conductor]
            info_text = (f"Area: {data['Aluminum_area_cmil']:,} cmil\n"
                        f"Stranding: {data['Stranding']}\n"
                        f"Diameter: {data['Outside_diameter_in']} in\n"
                        f"GMR: {data['GMR_ft']} ft")
            self.conductor_info.config(text=info_text)

    def receive_geometry_parameters(self, params):
        """Receive Dm from TLGD and update appropriate mode"""
        dm = params['dm']
        
        if self.mode_var.get() == "geometry":
            # Update Geometry Mode
            self.entry_Dm.config(state="normal")
            self.entry_Dm.delete(0, tk.END)
            self.entry_Dm.insert(0, f"{dm:.6f}")
            self.entry_Dm.config(state="normal")
            self.status_var.set(f"Dm = {dm:.3f}m received from geometry designer")
            
        elif self.mode_var.get() == "powerflow":
            # Update Power Flow Mode if using geometry
            if hasattr(self, 'geom_Dm'):
                self.geom_Dm.config(state="normal")
                self.geom_Dm.delete(0, tk.END)
                self.geom_Dm.insert(0, f"{dm:.6f}")
                self.geom_Dm.config(state="normal")
                self.status_var.set(f"Dm = {dm:.3f}m received from geometry designer")

    def calculate(self):
        try:
            self.result_text.delete(1.0, tk.END)
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            self.result_text.insert(tk.END, f"Analysis performed at: {timestamp}\n")
            self.result_text.insert(tk.END, "="*70 + "\n\n")
            
            if self.mode_var.get() == "geometry":
                self.calculate_geometry_mode()
            else:
                self.calculate_powerflow_mode()
                
            self.phasor_btn.config(state="normal")
                
        except Exception as e:
            messagebox.showerror("Calculation Error", f"An error occurred:\n{str(e)}")
            self.status_var.set(f"Error: {str(e)}")
            self.phasor_btn.config(state="disabled")

    def get_geometry_inputs(self):
        """Read and return geometry-related inputs in SI units (validated)."""
        conductor = self.conductor_combo.get()
        if not conductor or conductor not in aac_data:
            raise ValueError("Please select a valid conductor type")

        data = aac_data[conductor]
        try:
            n = int(self.entry_n.get())
        except Exception:
            raise ValueError("Invalid bundle count (n)")
        try:
            d = float(self.entry_spacing.get())
        except Exception:
            raise ValueError("Invalid bundle spacing")
        d_m = convert_length(d, self.spacing_unit.get(), "meters")
        try:
            Dm_m = float(self.entry_Dm.get())
        except Exception:
            # If entry is not set (designer not used), fallback to 5.0 m
            Dm_m = 5.0
        try:
            freq = float(self.entry_freq.get())
        except Exception:
            raise ValueError("Invalid frequency")
        try:
            length_km = float(self.entry_length.get())
        except Exception:
            raise ValueError("Invalid line length")
        param_units = self.param_units.get()
        output_unit = self.output_unit.get()

        return {
            "data": data,
            "n": n,
            "d_m": d_m,
            "Dm_m": Dm_m,
            "freq": freq,
            "length_km": length_km,
            "param_units": param_units,
            "output_unit": output_unit
        }

    def calculate_geometry_mode(self):
        """Calculate line parameters from geometry inputs"""
        # Read & validate inputs
        inputs = self.get_geometry_inputs()
        data = inputs["data"]
        n = inputs["n"]
        d_m = inputs["d_m"]
        Dm_m = inputs["Dm_m"]
        freq = inputs["freq"]
        length_km = inputs["length_km"]
        param_units = inputs["param_units"]
        output_unit = inputs["output_unit"]

        # Compute equivalent GMR for bundle (Ds)
        Ds = calc_Ds(n, data['GMR_ft'], d_m)

        # Compute R, XL, XC in requested param units
        R_output, XL_output, XC_output = calc_RLC_parameters(freq, Ds, Dm_m, data['Aluminum_area_cmil'], 20, param_units)

        # Also compute per-km values (used by ABCD)
        R_per_km, XL_per_km, XC_MΩ_km = calc_RLC_per_km(freq, Ds, Dm_m, data['Aluminum_area_cmil'])

        # Determine appropriate model and compute ABCD
        model_type, model_info = determine_line_model(length_km, freq)
        A, B, C, D = calc_ABCD_parameters(length_km, R_per_km, XL_per_km, XC_MΩ_km, freq, model_type)

        # Convert results for display
        Ds_out = convert_length(Ds, "meters", output_unit)
        Dm_out = convert_length(Dm_m, "meters", output_unit)

        unit_labels = {
            "per_km": ("Ω/km", "Ω/km", "MΩ·km"),
            "per_mile": ("Ω/mile", "Ω/mile", "MΩ·mile"),
            "per_m": ("Ω/m", "Ω/m", "MΩ·m"),
            "per_100km": ("Ω/100km", "Ω/100km", "MΩ·100km"),
            "per_100mile": ("Ω/100mile", "Ω/100mile", "MΩ·100mile")
        }
        r_unit, xl_unit, xc_unit = unit_labels.get(param_units, ("Ω/km", "Ω/km", "MΩ·km"))

        # Build result string
        result_str = f"""CONDUCTOR PARAMETERS:
  Conductor: {data['Aluminum_area_cmil']:,} cmil (selected)
  Bundle: {n} conductor(s)
  Equivalent GMR (Ds): {Ds_out:.8f} {output_unit}
  Geometric Mean Distance (Dm): {Dm_out:.6f} {output_unit}

LINE PARAMETERS ({param_units.replace('_', '-')}):
  Resistance (R): {R_output:.8f} {r_unit}
  Inductive Reactance (XL): {XL_output:.8f} {xl_unit}
  Capacitive Reactance (XC): {XC_output:.6f} {xc_unit}

LINE ANALYSIS:
  {model_info}
  Length: {length_km} km
  Frequency: {freq} Hz

ABCD PARAMETERS:
  A = {A:.8f}
  B = {B:.8f} Ω
  C = {C:.8f} S
  D = {D:.8f}
"""

        # Insert into results widget
        self.result_text.insert(tk.END, result_str)

        # Save last calculation state
        self.last_calculation = {
            'type': 'geometry',
            'Ds': Ds,
            'Dm': Dm_m,
            'R_per_km': R_per_km,
            'XL_per_km': XL_per_km,
            'XC_MΩ_km': XC_MΩ_km,
            'ABCD': (A, B, C, D),
            'length_km': length_km,
            'freq': freq,
            'Vs_phase': None,
            'Is_phase': None,
            'Vr_phase': None,
            'Ir_phase': None
        }

        self.status_var.set(f"Geometry analysis complete - {model_type} model used")
        self.result_text.insert(tk.END, "\n\n💡 Note: Phasor diagrams are available in Power Flow Mode")

    def calculate_powerflow_mode(self):
        """Calculate power flow parameters"""
        try:
            # Get line parameters based on input type
            if self.param_type.get() == "geometry":
                # Get geometry inputs for power flow
                conductor = self.geom_conductor_combo.get()
                if not conductor or conductor not in aac_data:
                    raise ValueError("Please select a valid conductor type")
                
                data = aac_data[conductor]
                n = int(self.geom_bundle.get())
                d_m = float(self.geom_spacing.get())
                Dm_m = float(self.geom_Dm.get())
                length_km = float(self.geom_length.get())
                freq = float(self.geom_freq.get())
                
                # Calculate parameters from geometry
                Ds = calc_Ds(n, data['GMR_ft'], d_m)
                R_per_km, XL_per_km, XC_MΩ_km = calc_RLC_per_km(freq, Ds, Dm_m, data['Aluminum_area_cmil'])
                
            elif self.param_type.get() == "RLC_per_unit":
                # Get per-unit parameters
                R_per_km = float(getattr(self, 'param_entry_0').get())
                XL_per_km = float(getattr(self, 'param_entry_1').get())
                XC_MΩ_km = float(getattr(self, 'param_entry_2').get())
                length_km = float(self.entry_line_length.get())
                freq = float(self.entry_freq_pf.get())
                
            elif self.param_type.get() == "RLC_total":
                # Get total parameters and convert to per-km
                R_total = float(getattr(self, 'param_entry_0').get())
                XL_total = float(getattr(self, 'param_entry_1').get())
                XC_total_MΩ = float(getattr(self, 'param_entry_2').get())
                length_km = float(self.entry_line_length.get())
                freq = float(self.entry_freq_pf.get())
                
                R_per_km = R_total / length_km
                XL_per_km = XL_total / length_km
                XC_MΩ_km = XC_total_MΩ / length_km
                
            else:  # YZ parameters
                Z_str = getattr(self, 'param_entry_0').get()
                Y_str = getattr(self, 'param_entry_1').get()
                length_km = float(self.entry_line_length.get())
                freq = float(self.entry_freq_pf.get())
                
                # Parse complex numbers
                Z_total = complex(Z_str.replace('j', 'j').replace('i', 'j'))
                Y_total = complex(Y_str.replace('j', 'j').replace('i', 'j'))
                
                # Convert to per-km
                Z_per_km = Z_total / length_km
                Y_per_km = Y_total / length_km
                
                R_per_km = Z_per_km.real
                XL_per_km = Z_per_km.imag
                # For capacitive reactance, we need to derive from susceptance
                XC_MΩ_km = 1 / (Y_per_km.imag * 1e6) if Y_per_km.imag != 0 else 1e6

            # Get power flow inputs
            Vr_kV = float(self.entry_voltage.get())
            Pr_MW = float(self.entry_power.get())
            pf = float(self.entry_pf.get())
            lead_lag = self.pf_var.get()
            
            # Determine model and calculate ABCD parameters
            model_type, model_info = determine_line_model(length_km, freq)
            A, B, C, D = calc_ABCD_parameters(length_km, R_per_km, XL_per_km, XC_MΩ_km, freq, model_type)
            
            # Calculate power flow
            results = calc_power_flow(Vr_kV, Pr_MW, pf, A, B, C, D, lead_lag)
            
            # Display results
            result_str = f"""POWER FLOW ANALYSIS:
{model_info}

LINE PARAMETERS:
  Length: {length_km} km
  Frequency: {freq} Hz
  Resistance: {R_per_km:.6f} Ω/km
  Inductive Reactance: {XL_per_km:.6f} Ω/km  
  Capacitive Reactance: {XC_MΩ_km:.6f} MΩ·km

RECEIVING END:
  Voltage: {Vr_kV} kV
  Power: {Pr_MW} MW
  Power Factor: {pf} ({'Lagging' if lead_lag == 'lag' else 'Leading'})

SENDING END RESULTS:
  Voltage: {results['Vs_kV']:.3f} kV
  Current: {results['Is_A']:.3f} A
  Power: {results['Ps_MW']:.3f} MW
  Reactive Power: {results['Qs_MVAR']:.3f} MVAR
  Efficiency: {results['efficiency']:.3f} %
  Voltage Regulation: {results['voltage_regulation']:.3f} %

ABCD PARAMETERS:
  A = {A:.8f}
  B = {B:.8f} Ω
  C = {C:.8f} S
  D = {D:.8f}
"""
            
            self.result_text.insert(tk.END, result_str)
            
            # Save for phasor diagram
            self.last_calculation = {
                'type': 'powerflow',
                'ABCD': (A, B, C, D),
                'length_km': length_km,
                'freq': freq,
                'Vs_phase': results['Vs_phase'],
                'Is_phase': results['Is_phase'],
                'Vr_phase': results['Vr_phase'],
                'Ir_phase': results['Ir_phase']
            }
            
            self.status_var.set("Power flow analysis complete")
            
        except Exception as e:
            raise ValueError(f"Power flow calculation error: {str(e)}")

    def plot_unified_phasors(self):
        """Plot unified phasor diagram using last calculation"""
        if self.last_calculation and self.last_calculation['type'] == 'powerflow':
            plot_unified_phasor_diagram(
                self.last_calculation['Vs_phase'],
                self.last_calculation['Is_phase'], 
                self.last_calculation['Vr_phase'],
                self.last_calculation['Ir_phase']
            )
        else:
            messagebox.showwarning("No Data", "Please run a power flow calculation first to generate phasor data.")

    def plot_transmission_geometry(self):
        """Open geometry designer in a new window"""
        geometry_window = tk.Toplevel(self.root)
        geometry_window.title("Transmission Line Geometry Designer")
        geometry_window.geometry("1200x800")
        
        designer = TransmissionLineGeometryDesigner(geometry_window, self)

    def clear_results(self):
        """Clear results text area"""
        self.result_text.delete(1.0, tk.END)
        self.status_var.set("Results cleared")

    def save_results(self):
        """Save results to file"""
        try:
            filename = f"transmission_line_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
            with open(filename, 'w') as f:
                f.write(self.result_text.get(1.0, tk.END))
            self.status_var.set(f"Results saved to {filename}")
            messagebox.showinfo("Save Successful", f"Results saved to {filename}")
        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save file: {str(e)}")

# ==============================
# === TEST SUITE AND MAIN ===
# ==============================

def run_automated_test_suite():
    """Run comprehensive test cases"""
    print("\n" + "="*80)
    print("🧪 TRANSMISSION LINE CALCULATOR - AUTOMATED TEST SUITE")
    print("="*80)
    
    test_results = []
    
    print("🚀 Running comprehensive test cases...")
    
    # Test Case 1: Basic Short Line
    print(f"\n1. Basic Short Line (50km, 230kV)")
    try:
        Ds = 0.0093 * 0.3048  # Waxwing
        Dm = 5.0
        R, XL, XC = calc_RLC_per_km(60, Ds, Dm, 266800)
        length = 50
        model_type, _ = determine_line_model(length, 60)
        A, B, C, D = calc_ABCD_parameters(length, R, XL, XC, 60, model_type)
        results = calc_power_flow(230, 100, 0.95, A, B, C, D, "lag")
        
        if 97 <= results['efficiency'] <= 99.5:
            test_results.append(("Basic Short Line", "PASS", f"η={results['efficiency']:.2f}%"))
            print(f"   ✅ PASS: Vs={results['Vs_kV']:.1f}kV, η={results['efficiency']:.2f}%")
        else:
            test_results.append(("Basic Short Line", "FAIL", f"η={results['efficiency']:.2f}%"))
    except Exception as e:
        test_results.append(("Basic Short Line", "FAIL", str(e)))

    # Test Case 2: Medium Line
    print(f"\n2. Medium Line (150km, 345kV)")
    try:
        Ds = 0.0108 * 0.3048  # Hen
        Dm = 7.0
        R, XL, XC = calc_RLC_per_km(60, Ds, Dm, 477000)
        length = 150
        model_type, _ = determine_line_model(length, 60)
        A, B, C, D = calc_ABCD_parameters(length, R, XL, XC, 60, model_type)
        results = calc_power_flow(345, 200, 0.92, A, B, C, D, "lag")
        
        if 96 <= results['efficiency'] <= 98.5:
            test_results.append(("Medium Line", "PASS", f"η={results['efficiency']:.2f}%"))
            print(f"   ✅ PASS: Vs={results['Vs_kV']:.1f}kV, η={results['efficiency']:.2f}%")
        else:
            test_results.append(("Medium Line", "FAIL", f"η={results['efficiency']:.2f}%"))
    except Exception as e:
        test_results.append(("Medium Line", "FAIL", str(e)))

    # Print Final Summary
    print("\n" + "="*80)
    print("📊 COMPREHENSIVE TEST SUMMARY")
    print("="*80)
    
    passed = sum(1 for result in test_results if result[1] == "PASS")
    failed = sum(1 for result in test_results if result[1] == "FAIL")
    total = len(test_results)
    
    print(f"Total Tests: {total}")
    print(f"✅ Passed: {passed}")
    print(f"❌ Failed: {failed}")
    print(f"📈 Success Rate: {passed/total*100:.1f}%")
    
    if failed > 0:
        print(f"\n🔍 FAILED TESTS:")
        for test_name, status, details in test_results:
            if status == "FAIL":
                print(f"   ❌ {test_name}: {details}")
    
    return test_results

if __name__ == "__main__":
    import sys
    
    # Run tests
    print("🚀 Initializing Transmission Line Calculator...")
    test_results = run_automated_test_suite()
    
    # Start GUI
    print("\n🎯 Starting GUI Application...")
    root = tk.Tk()
    app = TransmissionLineCalculator(root)
    root.mainloop()