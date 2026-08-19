import argparse
import os
import sys
import FreeCAD
import Import
import FreeCADGui

def render_step_to_png(step_path, png_path, width=1920, height=1080):
    # Initialize a headless document
    doc = FreeCAD.newDocument("RenderDoc")
    
    # Import the STEP file into the document
    print(f"Importing {step_path}...")
    Import.insert(step_path, doc.Name)
    
    # Recompute the document structure
    doc.recompute()
    
    # Initialize the GUI (necessary for rendering views)
    FreeCADGui.showMainWindow()
    from pivy import coin
    
    # Set the active document in GUI
    active_view = FreeCADGui.ActiveDocument.ActiveView
    
    # Fit all objects into the camera frame
    active_view.fitAll()
    
    # Optional: Set a standard isometric view orientation
    active_view.viewIsometric()
    
    # Render and save the image
    print(f"Saving render to {png_path}...")
    active_view.saveImage(png_path, width, height, "Current")
    print("Render complete!")

    # FreeCAD's Qt teardown crashes (or hangs indefinitely) after a headless
    # save; the image is on disk, so terminate without unwinding Qt.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)

# freecadcmd executes scripts with __name__ set to the file stem, not
# "__main__", so also run when launched by the FreeCAD console binary.
if __name__ == "__main__" or sys.argv[0].startswith("freecad"):
    parser = argparse.ArgumentParser(description="Render a STEP file to PNG")
    parser.add_argument("--step", default="model.step", help="input STEP file (default: model.step)")
    parser.add_argument("--png", default="output_render.png", help="output PNG file (default: output_render.png)")
    # freecadcmd leaves the raw command line in sys.argv: the script path plus
    # a --pass marker before every forwarded argument. Strip both so argparse
    # sees a clean argument list under freecadcmd and plain python alike.
    argv = [a for a in sys.argv[1:] if a != "--pass"]
    if argv and argv[0].endswith(".py"):
        argv = argv[1:]
    args = parser.parse_args(argv)

    render_step_to_png(args.step, args.png)
