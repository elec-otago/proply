import argparse
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
    
    # Clean up and close
    FreeCADGui.getMainWindow().close()
    print("Render complete!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render a STEP file to PNG")
    parser.add_argument("--step", default="model.step", help="input STEP file (default: model.step)")
    parser.add_argument("--png", default="output_render.png", help="output PNG file (default: output_render.png)")
    args = parser.parse_args()

    render_step_to_png(args.step, args.png)
