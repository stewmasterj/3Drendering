# view - Fortran fast 3D renderer to a linux framebuffer device

## features

Uses quaternion multiplication for smooth rotation of spatial vectors.
The scene is displayed using ray casting based on the spherical coordinates of the
 object's vertecies; rather than a projection onto a rectangular plane as most 3D
 engines do, this allows for a more realistic rendering for wide fields of view 
 without suffering from strange rectangular distorion, however it would appear as
 though looking through a fish-eye lense.  

The original invokation of this program is to specify an object to render and a
 configuration for the framebuffer device and display position.  
However, the current intention is to load objects from a `scene` file, and 
 specify that alone. Example scene files are found in the `scene` directory.
Example screenshots are in the `screenshot` directory. Example objects can be found
 in the `objFiles` directory. Programs that may help with object file formatting 
 and creation are located in `utils`.

### Requirements
may need to make yourself part of the 'video' group to write to '/dev/fb0'  
Depends on `fcurses`, `fbMod`, and `stringParseMods`. These can be found from the
 same github repository that this is from.


### Source files

-`renderMod.f90`  
	Uses `fbMod2.f90` for drawing to the screen and handles the objects and their positions in the scene.

-`view.f90`  
	Main routine for handling the orientation of the camera and any specific
 utility, this example has the capability for an asteroids game.

-`shareMods.f90`  
	This small module is used to share variables between `renderMod` and `view`.

## Internal Commandline Options

The `view` program allows for an internal commandline capability when ":" is pressed.
The subroutine is called from `renderMod` to read the line and interpret the duty.  
These commands can also be used in the `scene.rc` files to setup the initial scene.

Commands that take parameters show their values when executed sans parms  

### Command Parameters  Description  
 help               Displays this help message  
 q,quit,exit        terminates the program  
 fbpath F           path to framebuffer device, must have write permissions  
 dumpfile F         path for dump file for keyboard screenshots  
 tmpdir F           path to a temporary directory for scratch files  
 interactive        puts program into interactive rather than single render  
 contact N          calculates contact conditions for dynamic objects; N is the  
                     type of contact algorithm, 0=off, 1=elastic-Frictionless  
 FOV x y            field of view angles (radians)  
 subscreen x1 y1 x2 y2 pixel ranges, top-left and bottom right, (row column)  
 timestep r         time between frames for animation  
 dtheta r           angle change for key press  
 dspace r           distance increment for key press  
 width n            pixel width of screen  
 height n           pixel height of screen  
 timeout n          animation steps to run  
 linelength N       frame buffer line length Bytes  
 NobjectBuff N      allocate the object buffer to this many objects.  
                     WARNING will delete all existing objects  
 LoadObject FILE    will load object from file FILE  
 camera_pos x y z   location of camera  
 camera_orient a u v w  angle-axis orientation of camera  
 universe c x y z   set a universal shape (c=box|sphere), size parameters  
 periodic           sets the universal size to be periodic else will rebound  
 fbclear            clear the pixel buffer  
 redraw             writes the current pixel buffer to the frame buffer  
 pause              stops rendering animation  
 run                resumes rendering animation  
 picture            take a screen shot  
 record             toggle recording mode (video)  
 impulseControl     set input key mode to impulse  
 spatialControl     set input key mode to spatial  
 cpo I J (K)        copy object ID I to J (through K)  
 nearest,closest    reports the nearest vertex relative sph position and obj  
 underAim           reports the obj, tri, and node at center of screen  
 o I c              object ID 'I' definition 'c'  
    object command 'c' can be:  
     save s          save this object as file name 's'  
     name s          set name of object 'I'  
     mode s          set render mode of object 'I', (point|wire|solid)  
     mass r          set object's total mass (for collision)  
     radius r        set object's average radius (for contact)  
     offset x y z    set position offset of object 'I'  
     velocity u v w  set velocity of object 'I'  
     orient a u v w  set angle and rotation axis of object 'I'  
     spin a u v w    set angle rate and rotation axis of object 'I'  
     scale r         scale all local positions of object 'I' by multiple 'r'  
     add J           add a node to triangle J of object 'I', default J=underAIM  
     point J x y z   set local position of point 'J' in object 'I'  
     disp J u v w    displace local position of point 'J' in object 'I'  
     color J r g b   set color of point 'J' in object 'I' (J=0 for all points)  
     smooth J k      set smoothness of point J in object I (J=0 for all points)   
 orand c i j  v(:)   randomizes objects i-j property 'c' with ranges v(:)  



