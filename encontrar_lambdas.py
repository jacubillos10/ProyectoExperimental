#!/usr/bin/env/python
#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import csv

numEspectro=17;
nombreArchivo='min_espectro_sol'+str(numEspectro)+'.csv';
Datos=np.genfromtxt(nombreArchivo,delimiter=',', skip_header=1, usecols=(0,1,2,4));
pix=Datos[:,0];
Int=Datos[:,1];
posiciones=Datos[:,2];
elemento=np.loadtxt(nombreArchivo,delimiter=',',skiprows=1,usecols=3,dtype=str);
l_onda=Datos[:,3];
pix_reg=pix[elemento=="O2"];
l_reg=l_onda[elemento=="O2"];
coeff=np.polyfit(pix_reg,l_reg,1);
#pix_rem=pix[elemento==""];
#l_rem=coeff[0]*pix_rem+coeff[1];
pos_rem=posiciones[elemento==""];
# print(np.c_[pos_rem,l_rem]);
pix_Fe=pix[elemento=="Fe"];
l_Fe=l_onda[elemento=="Fe"];
pos_Fe=posiciones[elemento=="Fe"];
pix_o2=pix[elemento=="O2"];
l_o2=l_onda[elemento=="O2"];
pos_o2=posiciones[elemento=="O2"];
fig1=plt.figure();
ax1=fig1.add_subplot(111);
ax1.scatter(pix_o2,l_o2,marker='x',label="Oxígeno");
ax1.scatter(pix_Fe,l_Fe,marker='o',label="Hierro");
#for i in range(len(pix_o2)):
#	ax1.annotate(pos_o2[i],(pix_o2[i],l_o2[i]),fontsize="x-small");
#fin for
#for i in range(len(pix_Fe)):
#	ax1.annotate(pos_Fe[i],(pix_Fe[i],l_Fe[i]),fontsize="x-small");
##fin for 
ax1.set_xlabel("Posición en pixeles");
ax1.set_ylabel("Longitud de onda en nm");
ax1.legend(loc="upper left",shadow=True, fontsize="x-small");
plt.savefig("lineas_identificadas_sol"+str(numEspectro)+".png");
plt.close();

grado=2;
coeff_cuad=np.polyfit(pix_o2,l_o2,grado);
g=grado;
l_o2_calc=0;
l_Fe_calc=0;
for i in range(grado+1):
	l_o2_calc=l_o2_calc+coeff_cuad[i]*(pix_o2**g);
	l_Fe_calc=l_Fe_calc+coeff_cuad[i]*(pix_Fe**g);
	g=g-1;
#fin for 
print(coeff_cuad);
residual_o2=l_o2_calc-l_o2;
residual_Fe=l_Fe_calc-l_Fe;

fig2=plt.figure(figsize=(12,10));
ax2=fig2.add_subplot(111);
ax2.scatter(pix_o2,residual_o2,marker='x',label="Oxígeno");
ax2.scatter(pix_Fe,residual_Fe,marker='o',label="Hierro");
for i in range(len(pix_o2)):
	ax2.annotate(int(pos_o2[i]),(pix_o2[i],residual_o2[i]),fontsize='large');
#fin for 
for i in range(len(pix_Fe)):
	ax2.annotate(int(pos_Fe[i]),(pix_Fe[i],residual_Fe[i]),fontsize='large');
#fin for 
ax2.set_title("Gráfica de longitudes de onda residuales en función del pixel",fontsize='x-large');
ax2.set_xlabel("Posición en pixeles",fontsize='x-large');
ax2.set_ylabel("Diferencia de longitud de onda en nm",fontsize='x-large');
ax2.legend(loc="upper right",shadow=True, fontsize="large");
plt.savefig("residuales_sol"+str(numEspectro)+".png");
plt.close();

desv_o2=np.std(residual_o2)/(len(residual_o2)**0.5);
desv_Fe=np.std(residual_Fe)/(len(residual_Fe)**0.5);
incertidumbre=((desv_o2**2)+(desv_Fe**2))**0.5;
long_mediaFe=np.average(l_Fe);
print("El delta lambda en sol_",numEspectro," es de: ",np.average(residual_Fe)-np.average(residual_o2)," nm");
print("La incertidumbre es de: ",incertidumbre," nm");
print("La longitud de onda promedio es: ",long_mediaFe);
