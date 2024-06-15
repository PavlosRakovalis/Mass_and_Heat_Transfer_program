#Importing Libraries
import math
import sys
import matplotlib.pyplot as plt


#Stable Variables
kinematic_viscosity = 1.5*10**-5 #m^2/s
air_density = 1.187 #kg/m^3
air_specific_heat = 1000 #J/kgK
Sphere_diameter = 0.1 #m
sphere_density = 1000 #kg/m^3
sphere_specific_heat = 4180 #J/kgK
air_thermal_diffusivity = 2.208 * 10**(-5) #m^2/s 
air_thermal_conductivity = 0.02588 #W/mK

#eyresi epifaneias sfairas 
epiphaneia_sphairas = 4 * 3.14 * (Sphere_diameter/2)**2 #m^2
#eyresi mazas sfairas 
ogkos_sfairas = (4/3) * 3.14 * (Sphere_diameter/2)**3 #m^3
maza_sfairas = sphere_density * ogkos_sfairas #kg
#evresi piknotitas nerou stin eleytheri roi 
piknotita_nerou_eleytheris_rois = 0 #kg/m^3 

T_initial_matrix = [30, 40]  # Freestream air temperatures in Celsius
air_speed = 1  # Freestream air velocity in m/s

for a in range(2):
    #changing variables
    Tair = T_initial_matrix[a] + 273  #Kelvin
    enthalpia_eksatmisis_nerou = 2400000 #J/kg
    Dwater = 2.505 * 10**(-5) #m^2/s
    gas_constant_water = 461.52  #J/ kg*K                                        


    def water_vapor_pressure(T): #Pascal
        return 0.0001*100000*math.exp(13.766) * math.exp(-5132/T)  #Pascal

    Reynolds = air_speed * (Sphere_diameter) / kinematic_viscosity #adiastatos            



    #################evresi sinteletwn sinagwgis mazas kai thermotitas##################

    Prandtl = kinematic_viscosity/ air_thermal_diffusivity #adiastatos  

    Nusselt = 0.664 * Reynolds**(0.5) * Prandtl**(1/3) #adiastatos

    Sintelestis_thermikis_sinagwgis = Nusselt * air_thermal_conductivity / Sphere_diameter #adiastatos

    Schmidt = kinematic_viscosity / Dwater #adiastatos

    Sherwood = 0.664 * Reynolds**(0.5) * Schmidt**(1/3) #adiastatos

    sintelestis_sinagwgis_mazas = Sherwood * Dwater / Sphere_diameter #adiastatos



    #######################################Defining the Loop###########################



    steps = 200000
    Total_Time = 20000 #sec
    dt = Total_Time/steps #sec
    Tsphere = Tair #Kelvin

    T_sphere_matrix = []


    #Sintelestis_thermikis_sinagwgis = 10.416880407727518


    for i in range(steps):


        #ypologismos pyknotitas nerou stin epifaneia tis sfairas
        piknotita_nerou_epifaneias = water_vapor_pressure(Tsphere)/(gas_constant_water * Tsphere) #kg/m^3

        dT = dt *( Sintelestis_thermikis_sinagwgis * epiphaneia_sphairas *(Tair - Tsphere) - sintelestis_sinagwgis_mazas * epiphaneia_sphairas * (piknotita_nerou_epifaneias - piknotita_nerou_eleytheris_rois) * enthalpia_eksatmisis_nerou ) / (maza_sfairas*sphere_specific_heat) #Kelvin

    
        Tsphere = Tsphere + dT
        T_sphere_matrix.append(Tsphere-273)
        

        
    Minutes_matrix = [i*dt/60 for i in range(steps)]
       
    plt.plot(Minutes_matrix, T_sphere_matrix, label=f'Tair = {T_initial_matrix[a]}°C, U = {air_speed} m/s')
    plt.xlabel('Time (Minutes)')
    plt.ylabel('Tsphere(Celsius)')


    #################Lisi me hrisi tis analogias Chilton and Colburn##################
    Tguess = 15 + 273 #Kelvin
    Moriako_varos_nerou = 18 #gr/mol
    Moriako_varos_aera = 29 #gr/mol 
    Lewis = air_thermal_diffusivity / Dwater #adiastatos
    for i in range(100):
        T_chilton_coldburn = Tair - enthalpia_eksatmisis_nerou* Moriako_varos_nerou * (water_vapor_pressure(Tguess) - 0) /( (air_specific_heat * Lewis**(2/3) * Moriako_varos_aera * 101325  )) #Kelvin
        Tguess = Tguess - (Tguess - T_chilton_coldburn) / 2
        if abs(Tguess - T_chilton_coldburn) < 0.00001:
            print(f"T_chilton_coldburn = {T_chilton_coldburn - 273.15:.2f}°C for Tair = {T_initial_matrix[a]}°C, U = {air_speed} m/s")
            break

    # Add horizontal line at temperature = T_chilton_coldburn
    plt.axhline(y=T_chilton_coldburn - 273.15, color='red', linestyle='--', label=f'T_chilton_Colburn for Tair = {T_initial_matrix[a]}°C, U = {air_speed} m/s')

    # Annotate the horizontal line
    plt.annotate(f'T_chilton_Colburn for Tair = {T_initial_matrix[a]}°C, U = {air_speed} m/s', xy=(Minutes_matrix[-1], T_chilton_coldburn - 273.15), xytext=(Minutes_matrix[-1] - 100, T_chilton_coldburn - 273.15 + 1),
                 arrowprops=dict(facecolor='black', shrink=0.05))

# Add legend
plt.legend()

plt.show()  


for i in range(steps):
    if T_sphere_matrix[i] < 1.001*T_sphere_matrix[steps-1]:
        print(f"stability time at t = {i*dt/60:.2f} minutes for Tair = {T_initial_matrix[a]}°C, U = {air_speed} m/s")
        print(f'stable tempearture = {T_sphere_matrix[steps-1]:.2f}°C for Tair = {T_initial_matrix[a]}°C, U = {air_speed} m/s')
        break
#print(ogkos_sfairas)
