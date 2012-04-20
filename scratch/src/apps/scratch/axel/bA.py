'''
Created on 23.09.2011

@author: axel
'''
'''
Created on 03.05.2011

@author: axel
'''
from numpy import linspace, sort, arccos, zeros, sign, trapz, sin
from numpy.random import rand
from scipy.stats import norm
from math import pi, e

from numpy import array
from quaducom.resp_func.cb_short_fiber import CBShortFiber
from matplotlib import pyplot as p

#Dateneingabe:
E = 200000.
A = 0.0241
Le = 5
#Lg = Lc + Le
tau = 1.8
Df = 0.175
U = linspace( 0, Le, 40 )
#Ulim = ( Lg ** 2 - Lc ** 2 ) * tau / ( 2 * E * A )
#Ulim = ( U * E * A * tau * pi * Df ) ** 0.5 / tau - Le

Cbsf = CBShortFiber()



#Heavyside
def H( x ):
    return sign( sign( x ) + 1 )

def pullout( U, phi, tau, f ):
    y = ( ( H( Ulim - U ) * ( U * E * A * tau * pi * Df ) ** 0.5 ) + H( U - Ulim ) * ( tau * ( Le - U ) ) ) * e ** ( f * phi )
    return y

#Tau Numerisch
qTAU = rand( 1000 )
ppfTAU = ( ( norm( 1.7677, 0.6195 ).ppf( qTAU ) ) ** 2 ) ** 0.5
pdfTAU = norm( 1.7677, 0.6195 ).pdf( ppfTAU )

#Phi Numerisch
qPHI = rand( 1000 )
ppfPHI = arccos( -2 * qPHI + 1 ) / 2
pdfPHI = sin( 2 * ppfPHI )

#Beta Numerisch
#qBETA = rand( 70 )
#ppfBETA = ( ( norm( 0.0635, 0.0448 ).ppf( qBETA ) ) ** 2 ) ** 0.5
#pdfBETA = norm( 0.0635, 0.0448 ).pdf( ppfBETA )

#Pmean,Pvar
list_PmeanTAU_num = []
list_Pvar_num = []
list_PmeanMINUSPsdt = []
list_PmeanPLUSPsdt = []

#PHI
list_PmultiPHI_num = []
list_PmeanPHI_num = []

#Beta
list_PmultiBETA_num = []
#list_PmeanBETA_num = []
list_PvarPHI_num = []
list_PmultivarPHI_num = []
#list_PmultivarBETA_num = []
#list_PvarBETA_num = []

for Ul in U:
    #for BETA in sort( ppfBETA ):
        for PHI in sort( ppfPHI ):
            #Mean TAU
            ymeanTAU = Cbsf( Ul, sort( ppfTAU ), Le, Df, E, 5, PHI, 0.87, 0, 0, 1e15 ) * sort( pdfTAU )
            PmeanTAU_num = trapz( ymeanTAU, sort( ppfTAU ) )
            ymeanPHI = PmeanTAU_num * sin( PHI )
            list_PmultiPHI_num.append( ymeanPHI )
            #Var TAU
            yvarTAU = ( Cbsf( Ul, qTAU, Le, Df, E, 5, PHI, 0.87, 0, 0, 1e15 ) - ymeanTAU ) ** 2 * sort( pdfTAU )
            PvarTAU_num = trapz( yvarTAU, sort( ppfTAU ) )
            yvarPHI = PvarTAU_num * sin( PHI )
            list_PmultivarPHI_num.append( yvarPHI )
        #MEAN PHI
        PmeanPHI_num = trapz( list_PmultiPHI_num, sort( ppfPHI ) )
        list_PmultiPHI_num = []
        list_PmeanPHI_num.append( PmeanPHI_num )
        #VAR PHI
        PvarPHI_num = trapz( list_PmultivarPHI_num, sort( ppfPHI ) )
        list_PmultivarPHI_num = []
        list_PvarPHI_num.append( PvarPHI_num )
        #MEAN Beta
        #ymeanBETA = PmeanTAU_num * norm( 0.0635, 0.0448 ).pdf( BETA )
        #list_PmultiBETA_num.append( ymeanBETA )
        #Var Beta
        #yvarBETA = PvarTAU_num * norm( 0.0635, 0.0448 ).pdf( BETA )
        #list_PmultivarBETA_num.append( yvarBETA )
    #MEAN
    #PmeanBETA_num = trapz( list_PmultiBETA_num, sort( ppfBETA ) )
    #list_PmultiBETA_num = []
    #list_PmeanBETA_num.append( PmeanBETA_num )
    #VAR
    #PvarBETA_num = trapz( list_PmultivarBETA_num, sort( ppfBETA ) )
    #list_PmultivarBETA_num = []
    #list_PvarBETA_num.append( PvarBETA_num )




Pmean = array( list_PmeanPHI_num )
Pvar = array( list_PvarPHI_num )
print Pmean, 0, 0, 0, 0, 0, Pvar

Array_PplusStd = Pmean + Pvar ** 0.5
Array_PminusStd = Pmean - Pvar ** 0.5

#print len( list_PmeanBETA_num ), len( list_PvarBETA_num )
#print list_Pvar_num
x = U
y1 = list_PmeanPHI_num
y2 = Array_PplusStd
y3 = Array_PminusStd
p.plot( x, y1, x, y2, 'k--', x, y3, 'k--' )
p.title( 'numeric determining of mean pull out and scatterband' )
p.xlabel( '[mm]' )
p.ylabel( '[N]' )
#x,y2,'k--',x,y3,'k--'

#y=pullout(5.0,0,x_PDF_TAU,0.0635,0.8706)*pdfTAU
#Pmean_num=trapz(y,x_PDF_TAU)
#print Pmean_num
#p.plot(x_PDF_TAU,pdfTAU)
#p.plot(x,y)
#p.show()
'''
P_array=zeros(len(U))
for TAU in ppfTAU:   
        y=pullout(U,0,TAU,0.0635,0.8706)
        x=U
        P_array+=y    
Pmean_monte=P_array/len(ppfTAU)
print len(Pmean_monte)
print len(ppfTAU)

Pvar_array=zeros(len(U))
for TAU in ppfTAU:
    y=(pullout(U,0,TAU,0.0635,0.8706)-Pmean_monte)**2
    Pvar_array+=y
Pvar_monte=Pvar_array/len(ppfTAU)
print len(Pvar_monte)
Punder_monte=Pmean_monte-(Pvar_monte**0.5)
Pupper_monte=Pmean_monte+(Pvar_monte**0.5)
x=U
y4=Pmean_monte
y5=Punder_monte
y6=Pupper_monte
p.plot(x,y4,x,y5,x,y6)
'''
p.show()
