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
    # Define your file paths here or pass them as arguments
    input_step = "model.step"
    output_png = "output_render.png"
    
    render_step_to_png(input_step, output_png)
