#declare CameraBlur = 1; //turns camera blur on and off, 1 = on
#declare AreaLight = 1; //turns softening of shadows on and off, 1 = on

#include "colors.inc"

#declare lookat = <9.5,1,9.5>;
#declare focalpoint = <9.5,1,9.5>;

//focal blur camera
camera {
  location 50*<-.5, .99,  -1>
  look_at   lookat
  angle 20
  right     x*image_width/image_height
  #if (CameraBlur)
  	aperture 3          // [0...N] larger is narrower depth of field (blurrier)
  	blur_samples 1000        // number of rays per pixel for sampling
  	focal_point focalpoint    // point that is in focus <X,Y,Z>
  	confidence 0.95           // [0...<1] when to move on while sampling (smaller is less accurate)
  	variance 1/200            // [0...1] how precise to calculate (smaller is more accurate)
  #end
}

//light source
light_source {
  0*x                 // light's position (translated below)
  color rgb 1.0       // light's color
  #if (AreaLight)
  	area_light
  	<200, 0, 0> <0, 200, 0> // lights spread out across this distance (x * z)
  	4, 4                // total number of lights in grid (4x*4z = 16 lights)
  	adaptive 3          // 0,1,2,3...
  	jitter              // adds random softening of light
  	circular            // make the shape of the light circular
  	orient              // orient light
  #end
  translate 200*<-.2, .7,  -1>
}


// set a color of the background (sky)
background { color White }

#declare FOVrad = .03;

#declare FOVplane = union{
	sphere{<.5,1,.5>,FOVrad}
	sphere{<16.5,1,.5>,FOVrad}
	sphere{<16.5,1,16.5>,FOVrad}
	sphere{<.5,1,16.5>,FOVrad}
	cylinder{<.5,1,.5>,<16.5,1,.5>,FOVrad}
	cylinder{<16.5,1,.5>,<16.5,1,16.5>,FOVrad}
	cylinder{<16.5,1,16.5>,<.5,1,16.5>,FOVrad}
	cylinder{<.5,1,16.5>,<.5,1,.5>,FOVrad}
}

object{FOVplane pigment{ color Black} no_shadow no_reflection}

#declare F_ODF = finish {
    	ambient .4
    	specular .4
  	}
  	
#declare T_ODF = texture {
 	finish{F_ODF}
    	pigment{rgb <1,1,1>}
  	}

#declare ODF = mesh2 {
  vertex_vectors {
		#fopen XYZdat "ODF_verts.txt" read
		#read (XYZdat,DatLength)
		DatLength,
		#declare Count = 0;	
		#while (Count<DatLength-1)
			#read (XYZdat,Xval,Yval,Zval)
			<Xval, Zval, Yval>,
			#declare Count = Count+1;	
		#end
		#read (XYZdat,Xval,Yval,Zval)
		<Xval, Zval, Yval>
		#fclose XYZdat
  }
  
  	normal_vectors {
		#fopen XYZdat "ODF_norms.txt" read
		#read (XYZdat,DatLength)
		DatLength,
		#declare Count = 0;	
		#while (Count<DatLength-1)
		#read (XYZdat,Xval,Yval,Zval)
			<Xval, Zval, Yval>,
			#declare Count = Count+1;	
		#end
		#read (XYZdat,Xval,Yval,Zval)
		<Xval, Zval, Yval>
		#fclose XYZdat
	}
	  
	  texture_list {
		#fopen XYZdat "ODF_c.txt" read
		#read (XYZdat,DatLength)
		DatLength,
		#declare Count = 0;	
		#while (Count<DatLength-1)
		#read (XYZdat,Xval,Yval,Zval)
			texture{finish{F_ODF} pigment{rgb <Xval, Yval, Zval>}}
			#declare Count = Count+1;	
		#end
		#read (XYZdat,Xval,Yval,Zval)
		texture{finish{F_ODF} pigment{rgb <Xval, Yval, Zval>}}
		#fclose XYZdat
	   }

	face_indices {
		#fopen XYZdat "ODF_tri.txt" read
		#read (XYZdat,DatLength)
		DatLength,	
		#declare Count = 0;	
		#while (Count<DatLength-1)
			#read (XYZdat,Xval,Yval,Zval)
			<Xval, Yval, Zval>, Xval, Yval, Zval
			#declare Count = Count+1;
		#end
		#read (XYZdat,Xval,Yval,Zval)
		<Xval, Yval, Zval>, Xval, Yval, Zval
		#fclose XYZdat
	}  
}

object{ODF scale <1,1,1> rotate 0*z rotate -0*y no_shadow no_reflection}



