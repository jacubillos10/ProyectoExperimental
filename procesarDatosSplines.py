#!/usr/bin/env/python
#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import optimize

numEspectro=17;
nPuntos_ID=50;
derivada2_minima=0;
nombreArchivo='sol'+str(numEspectro)+'.dat';
Datos=np.genfromtxt(nombreArchivo,delimiter=' ');
x=Datos[:,0];
I=Datos[:,1];

def splin(x,datos_x,datos_y):
	p=np.zeros(len(datos_x));
	for i in range(len(datos_x)):
		num=1;
		den=1;
		for j in range(len(datos_x)):
			if i==j:
				pass;
			else:
				num=num*(x-datos_x[j]);
				den=den*(datos_x[i]-datos_x[j]);
			#fin if
		#fin for
		p[i]=datos_y[i]*(num/den);
	#fin for
	Suma=sum(p);
	return Suma;
#fin derivada_poli

#Genera la matriz que nos va a permitir sacar las derivadas
#por favor, coloque un número par. Esta función generará una matriz nxn
#con un n dado por el usuario
def generar_matriz(n):
	Matriz_raw=np.zeros((n,n));
	for j in range(n):
		if j<(n/2):
			indice_k=j-(n/2);
		else:
			indice_k=j-(n/2)+1;
		#fin if 
		for l in range(n):
			Matriz_raw[j,l]=(indice_k**(l+1))/np.math.factorial(l+1);
		#fin for 
	#fin for 
	Matriz=np.linalg.inv(Matriz_raw);
	return Matriz;
#fin función

#halla la m derivada de cada uno de los datos de entrada
#Por favor, coloque un m<=n, por ejemplo si m=2 halla la 
#segunda derivada de los datos usando un orden n de Taylor
def generar_datos_derivada_m(m,n,x,y):
	resp=np.zeros(len(x));
	Matriz_der=generar_matriz(n);
	for i in range(len(x)):
		atras=len(x)-1-i;
		if (i<(n/2)) or (atras<(n/2)):
			resp[i]=0;
		else:
			v_ref=np.zeros(n);
			for j in range(n):
				if j<(n/2):
					indice_k=j-int(n/2);
				else:
					indice_k=j-int(n/2)+1;
				#fin if
				v_ref[j]=y[i+indice_k]-y[i];
			#fin for 
			vec_mult=Matriz_der[m-1,:];
			prod=vec_mult.dot(v_ref);
			h=(x[i+1]-x[i])**(m-1);
			resp[i]=prod/h;
		#fin if
	#fin for 
	return resp;
#fin función

def encontrar_menores_intensidades(x,I,tam,vecinos):
	I_menores=max(I)*np.ones(tam);
	indices=np.zeros(tam);
	xs=np.zeros(tam);
	Is=np.zeros(tam);
	for k in range(len(x)):
		termino=False;
		l=0;
		while l<tam and termino==False:
			if I[k]<I_menores[l]:
				if abs(x[k]-x[int(indices[l])])<=vecinos:
					I_menores[l]=I[k];
					indices[l]=k;
					termino=True;
				else:
					if l<(tam-1):
						for g in range(tam-1-l):
							I_menores[-g-1]=I_menores[-g-2];
							indices[-g-1]=indices[-g-2];
						#fin for 
					#fin if 
					I_menores[l]=I[k];
					indices[l]=k;
					termino=True;
				#fin if 
			else:
				if abs(x[k]-x[int(indices[l])])<=vecinos:
					termino=True;
				#fin if 
			#fin if 
			l=l+1;
		#fin while
	#fin for 
	for n in range(tam):
		xs[n]=x[int(indices[n])];
		Is[n]=I[int(indices[n])];
	#fin for
	posiciones=np.argsort(xs);
	xs1=np.zeros(len(xs));
	Is1=np.zeros(len(xs));
	for k in range(len(xs)):
		xs1[k]=xs[posiciones[k]];
		Is1[k]=Is[posiciones[k]];
	#fin for
	return [xs1, Is1];
#fin función

def filtro1(x,I,der2,der_min):
	x_fin=[];
	I_fin=[];
	for k in range(len(x)):
		if k>0 and k<(len(x)-1):
			if I[k-1]>I[k] and I[k+1]>I[k] and der2[k]>der_min:
				x_fin.append(x[k]);
				I_fin.append(I[k]);
			#fin if 
		#fin if
	#fin for
	x_fin=np.array(x_fin);
	I_fin=np.array(I_fin);
	return [x_fin,I_fin];
#fin filtro 1

der_2=generar_datos_derivada_m(2,16,x,I);
der_min=derivada2_minima;
Filtrados1=filtro1(x,I,der_2,der_min);
x_filtro=Filtrados1[0];
I_filtro=Filtrados1[1];
xs_Is=encontrar_menores_intensidades(x_filtro,I_filtro,nPuntos_ID,3);

x_fin=xs_Is[0];
I_fin=xs_Is[1];
pixeles=np.zeros(len(x_fin));
Iopt=np.zeros(len(x_fin));
ajustes_x=[];
ajustes_I=[];
puntos_x=[];
puntos_I=[];
for j in range(len(x_fin)):
	maxVec=6;
	vec=1;
	act=0;
	termino=False;
	pos=int(x_fin[j]);
	while vec<maxVec and termino==False:
		if I[pos+vec]>I[pos+act] and I[pos-vec]>I[pos-act]:
			vec=vec+1;
			act=act+1;
		else:
			termino=True;
		#fin if
	#fin while
	xs=np.linspace(pos-vec,pos+vec,2*vec+1, dtype=int);
	Is=np.zeros(2*vec+1);
	for l in range(2*vec+1):
		Is[l]=I[xs[l]];
	#fin for 
	Resultados=np.polyfit(xs,Is,4); 
	Finura=100;
	xsp=np.linspace(min(xs)+1,max(xs)-1,Finura);
	Isp=np.zeros(len(xsp));
	funcSplin=lambda x: splin(x,xs,Is);
	opti=optimize.minimize(funcSplin,pos);
	pixel_optimo=opti.x[0];
	pixeles[j]=round(pixel_optimo,4);
	Iopt[j]=round(opti.fun,4);
	for l in range(len(xsp)):
		Isp[l]=splin(xsp[l],xs,Is);
	#fin for 
	puntos_x.append(xs);
	puntos_I.append(Is);
	ajustes_x.append(xsp);
	ajustes_I.append(Isp);
#fin for 

 
fig1=plt.figure(figsize=(14,9));
ax1=fig1.add_subplot(111);
ax1.plot(x,I);
ax1.scatter(pixeles,Iopt,color='orange');
for n in range(len(ajustes_x)):
	#ax1.scatter(puntos_x[n],puntos_I[n]);
	ax1.plot(ajustes_x[n], ajustes_I[n], color='purple');
#fin for
pos1=np.argsort(pixeles);
pixeles_ord=np.zeros(len(pixeles));
Int_ord=np.zeros(len(pixeles));
nums=np.zeros(len(pixeles));
#lamb_tent=np.zeros(len(pixeles));
for i in range(len(pixeles)):
	nums[i]=i;
	pixeles_ord[i]=pixeles[pos1[i]];
	Int_ord[i]=Iopt[pos1[i]];
	#lamb_tent[i]=round(((633-626)/1530)*pixeles_ord[i]+626,2);
	ax1.annotate(i,(pixeles_ord[i],Int_ord[i]));
#fin for 
#ax1.set_xlim([68,83]);
#ax1.set_ylim([18000,32000]);
ax1.set_xlabel("Numero de pixel en 'x'", fontsize='xx-large');
ax1.set_ylabel("Intensidad (unidades arbitrarias)", fontsize='xx-large');
ax1.set_title("Intensidad en función de la posición del pixel", fontsize='xx-large');
plt.savefig("grafica_del_espectro_sol"+str(numEspectro)+".png");

print("Los pixeles con I mínimo medido es: ", x_fin);
print("Los pixeles con I mínimo ajustado es: ", pixeles);
print("Los valores de intensidad son: ",Iopt);
encabezado=['posicion pixel', 'I minima','numero_pixel'];
datos_exportacion=np.c_[pixeles_ord,Int_ord,nums];
nombreCSV='min_espectro_sol'+str(numEspectro)+'.csv';
with open(nombreCSV, 'w', encoding='UTF8', newline='') as f:
	writer=csv.writer(f);
	writer.writerow(encabezado);
	writer.writerows(datos_exportacion);
#fin with 
