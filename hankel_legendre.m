(* ::Package:: *)

(* ::Section:: *)
(*Hankel function mapped on to the surface of a sphere*)


(* ::Text:: *)
(*Spherical Harmonics and the two components of the Hankel function*)


(* ::Input:: *)
(*l=20;*)
(*m=5;*)
(*f=Re;*)
(*gradients=ColorData["RedGreenSplit"];*)
(*SphericalPlot3D[1,{\[Theta],0,\[Pi]},{\[Phi],0,2 \[Pi]},*)
(*ColorFunction->(gradients[.5+f[SphericalHarmonicY[l,m,#4,#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*SphericalPlot3D[1,{\[Theta],0,\[Pi]},{\[Phi],0,2 \[Pi]},*)
(*ColorFunction->(gradients[.5+f[(HankelH1[m,Sqrt[l(l+1)]*#4] )*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*SphericalPlot3D[1,{\[Theta],0,\[Pi]},{\[Phi],0,2 \[Pi]},*)
(*ColorFunction->(gradients[.5+f[(HankelH2[m,Sqrt[l(l+1)]*#4])*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)


(* ::Section:: *)
(*Legendre function*)


(* ::Text:: *)
(*Spherical Harmonics and the two components of the Legendre function*)


(* ::Input:: *)
(*l=20;*)
(*m=5;*)
(*f=Re;*)
(*gradients=ColorData["RedGreenSplit"];*)
(*SphericalPlot3D[1,{\[Theta],0,\[Pi]},{\[Phi],0,2 \[Pi]},*)
(*ColorFunction->(gradients[.5+f[SphericalHarmonicY[l,m,#4,#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*SphericalPlot3D[1,{\[Theta],0,\[Pi]},{\[Phi],0,2 \[Pi]},*)
(*ColorFunction->(gradients[.5+f[Sqrt[(2l+1)/(2Pi) (l-m)!/(l+m)!]*(LegendreP[l,m,Cos[#4 ]]+ ((2*I)/Pi)LegendreQ[l,m,Cos[#4] ])*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*SphericalPlot3D[1,{\[Theta],0,\[Pi]},{\[Phi],0,2 \[Pi]},*)
(*ColorFunction->(gradients[.5+f[Sqrt[(2l+1)/(2Pi) (l-m)!/(l+m)!]*(LegendreP[l,m,Cos[#4 ]]- ((2*I)/Pi)LegendreQ[l,m,Cos[#4] ])*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)



