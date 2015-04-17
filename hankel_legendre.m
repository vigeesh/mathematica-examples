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


(* ::Section:: *)
(*Selecting a small patch for the analysis*)


(* ::Text:: *)
(*Selecting small patch on a sphere for analysis*)


(* ::Input:: *)
(*l=20;*)
(*m=5;*)
(*f=Re;*)
(*\[Theta]c=90;*)
(*d\[Theta] = 20;*)
(*\[Phi]c=0;*)
(*d\[Phi] = 100;*)
(*\[Theta]c=\[Theta]c*\[Pi]/180;*)
(*d\[Theta] =d\[Theta] *\[Pi]/180;*)
(*\[Phi]c=\[Phi]c*\[Pi]/180;*)
(*d\[Phi] = d\[Phi]*\[Pi]/180;*)
(*gradients=ColorData["RedGreenSplit"];*)
(*(*The Spherical Harmonics*)*)
(*SphericalPlot3D[1,{\[Theta],\[Theta]c-d\[Theta]/2,\[Theta]c+d\[Theta]/2},{\[Phi],\[Phi]c-d\[Phi]/2,\[Phi]c+d\[Phi]/2},*)
(*ColorFunction->(gradients[.5+f[SphericalHarmonicY[l,m,#4,#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*(*Hankel A*)*)
(*SphericalPlot3D[1,{\[Theta],\[Theta]c-d\[Theta]/2,\[Theta]c+d\[Theta]/2},{\[Phi],\[Phi]c-d\[Phi]/2,\[Phi]c+d\[Phi]/2},*)
(*ColorFunction->(gradients[.5+f[(HankelH1[m,Sqrt[l(l+1)]*#4] )*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*(*Hankel A*)*)
(*SphericalPlot3D[1,{\[Theta],\[Theta]c-d\[Theta]/2,\[Theta]c+d\[Theta]/2},{\[Phi],\[Phi]c-d\[Phi]/2,\[Phi]c+d\[Phi]/2},*)
(*ColorFunction->(gradients[.5+f[(HankelH2[m,Sqrt[l(l+1)]*#4])*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*(*Legendre A*)*)
(*SphericalPlot3D[1,{\[Theta],\[Theta]c-d\[Theta]/2,\[Theta]c+d\[Theta]/2},{\[Phi],\[Phi]c-d\[Phi]/2,\[Phi]c+d\[Phi]/2},*)
(*ColorFunction->(gradients[.5+f[Sqrt[(2l+1)/(2Pi) (l-m)!/(l+m)!]*(LegendreP[l,m,Cos[#4 ]]+ ((2*I)/Pi)LegendreQ[l,m,Cos[#4] ])*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
(*(*Legendre B*)*)
(*SphericalPlot3D[1,{\[Theta],\[Theta]c-d\[Theta]/2,\[Theta]c+d\[Theta]/2},{\[Phi],\[Phi]c-d\[Phi]/2,\[Phi]c+d\[Phi]/2},*)
(*ColorFunction->(gradients[.5+f[Sqrt[(2l+1)/(2Pi) (l-m)!/(l+m)!]*(LegendreP[l,m,Cos[#4 ]]- ((2*I)/Pi)LegendreQ[l,m,Cos[#4] ])*Exp[I*m*#5]]]&),Mesh->{17},MeshStyle->Opacity[0.2],*)
(*Boxed->False,Axes->False,ColorFunctionScaling->False,PlotPoints->200,SphericalRegion->True,ViewAngle->.3,ImageSize->100]*)
