README
Computer Graphics Raytrace Renderer

The main program is located in /assignment04/src/apps/04_pathtrace.cpp

pathtrace.cpp was built on a series of preceeding labs as well as framework provided by the instructor
including image.cpp, vmath.h, and scene.h


Project Build Instructions from command line:
	Navigate to /assignment04
	>mkDir Xcode
	>cd Xcode
	>cmake -G "Xcode" ..

Running the Project:
	Open /assignment04/Xcode/04pathtrace.xcodeproj using Xcode
	Provide the required json parameter by editing Product->Scheme->Edit Scheme
		In "Argument Passed on Launch" under the "Arguments" tab add the full path name to a json file located in /assignment04/scenes/
		ex.
			/Users/dackerman/Desktop/portfolio/clemson/code/Ackerman_Raytrace/assignment04/scenes/07_cb_indirect.json

		Under the "Info" tab in the Edit Scheme menu choose 04_pathtrace as the executable file.

	Run program using Product->Run. The rendered image will save as a png in /assignment04/scenes/
	The name of the rendered image will be the same as the json argument
		ex.
			07_cb_indirect.png



