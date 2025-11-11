# Laplacian Deformation GUI - User Guide

## Overview

The GUI application provides an interactive interface for performing Laplacian mesh deformation on 3D models. You can select fixed and movable handles, drag them to desired positions, and see the mesh deform smoothly while preserving local geometric details.

## How to Run

```powershell
cd build\bin
.\minimeshgui.exe ..\..\meshes\woody-lo.obj
```

Or for any mesh file:
```powershell
.\minimeshgui.exe path\to\your\mesh.obj
```

## Interface Overview

The application consists of two windows:
1. **3D Viewer Window** - Displays the mesh and allows interaction
2. **Controls Window** - Contains all buttons and settings

## Step-by-Step Usage

### Step 1: Enable Deformation Mode

1. Click **"Toggle Deformation Mode"** button in the Controls window
2. You should see a message in the console confirming deformation mode is enabled
3. The mouse mode automatically switches to "Select vertex"

### Step 2: Add Fixed Handles (Red)

Fixed handles are vertices that stay in place during deformation.

1. **Click** on vertices you want to fix
2. Selected vertices will appear with **RED spheres**
3. Typical use: Click vertices at the bottom or base of the model
4. Click the same vertex again to remove it as a handle

### Step 3: Add Movable Handles (Green)

Movable handles are vertices you will drag to deform the mesh.

1. **Hold Shift + Click** on vertices you want to move
2. Selected vertices will appear with **GREEN spheres**
3. Typical use: Click vertices at the top or areas you want to deform
4. Shift+Click the same vertex again to remove it

### Step 4: Choose Weight Type

In the "Laplacian Deformation" panel:
- **Cotangent weights** (default) - More accurate, geometry-aware
- **Uniform weights** - Simpler, faster, but less accurate

### Step 5: Move Handles

1. Switch mouse mode to **"Move vertex"** in the Controls window
2. Click and drag any **GREEN handle** (movable vertex)
3. The vertex will move interactively as you drag
4. You can drag multiple times to accumulate displacement

### Step 6: Apply Deformation

1. Click **"Apply Deformation"** button
2. The algorithm will:
   - Build the Laplacian matrix
   - Compute Laplacian coordinates
   - Apply constraints
   - Solve the linear system
   - Update the mesh

3. Watch the console for progress and diagnostics:
   ```
   === Applying Laplacian Deformation ===
   Fixed handles: 5
   Movable handles: 3
   Weight type: cotangent
   Building Laplacian matrix...
   ...
   === Diagnostics ===
   Constraint error (RMSE): 0.00e+00
   Residual norm (RMS): 1.23e-13
   ```

### Step 7: Iterate or Reset

After applying deformation, you can:
- **Add more handles** and apply again for further deformation
- **Clear Handles** to start over with handle selection
- **Reset Mesh** to restore the original shape

## Controls Panel Reference

### View Mesh Components
- ☑ **Show vertices** - Display vertex points
- ☑ **Show edges** - Display mesh edges
- ☑ **Show faces** - Display mesh faces (default on)
- ☑ **Show axis** - Display coordinate axes
- ☑ **Show lighting** - Enable/disable lighting

### Mouse Functionality
- **Pan and zoom** - Navigate the 3D view (rotate, zoom, pan)
- **Select vertex** - Click to select handles (use in deformation mode)
- **Move vertex** - Drag movable handles to displace them

### Laplacian Deformation Panel
- **Toggle Deformation Mode** - Enable/disable handle selection
- **Cotangent weights / Uniform weights** - Choose weight type
- **Apply Deformation** - Execute the deformation algorithm
- **Clear Handles** - Remove all fixed and movable handles
- **Reset Mesh** - Restore original mesh geometry

## Tips and Best Practices

### Handle Placement

1. **Fixed Handles (Red)**:
   - Place at boundaries or areas you want to anchor
   - Bottom vertices for objects standing upright
   - Base of limbs or appendages
   - Minimum 1-2 fixed handles required

2. **Movable Handles (Green)**:
   - Place at areas you want to deform
   - Fewer handles = smoother deformation
   - More handles = more control, less freedom
   - Try 1-5 movable handles for best results

### Deformation Strategy

1. **Start Simple**: Begin with just a few handles
2. **Incremental Changes**: Apply small deformations, then add more handles
3. **Anchor First**: Set fixed handles before movable ones
4. **Test and Iterate**: Apply, observe, adjust, repeat

### Common Workflows

**Bending**: 
- Fix bottom vertices (red)
- Move top vertices sideways (green)

**Stretching**:
- Fix bottom vertices (red)
- Move top vertices upward (green)

**Twisting**:
- Fix bottom center (red)
- Move top vertices in circular pattern (green)

**Local Editing**:
- Fix surrounding vertices (red)
- Move central area (green)

## Keyboard Shortcuts

- **Mouse Wheel** - Zoom in/out
- **Left Drag** (Pan mode) - Rotate view
- **Right Drag** (Pan mode) - Pan view
- **Shift** - Hold while clicking to add movable handles

## Diagnostics Interpretation

After applying deformation, check the console output:

- **Constraint error**: Should be ~0 (constraints satisfied exactly)
- **Residual norm**: Should be < 1e-10 (good numerical accuracy)

If you see:
- High constraint error (> 1e-6): Something went wrong
- High residual (> 1e-5): Numerical issues, try fewer handles

## Troubleshooting

**Problem**: Can't select vertices
- **Solution**: Make sure "Toggle Deformation Mode" is ON and mouse mode is "Select vertex"

**Problem**: Vertices don't move when dragging
- **Solution**: Switch mouse mode to "Move vertex" and only drag GREEN handles

**Problem**: Deformation looks wrong or mesh explodes
- **Solution**: 
  - Use "Reset Mesh" to start over
  - Try fewer movable handles
  - Use cotangent weights instead of uniform
  - Make sure you have at least one fixed handle

**Problem**: Apply Deformation button does nothing
- **Solution**: You need at least one handle (fixed or movable)

**Problem**: Spheres (handles) not visible
- **Solution**: Make sure "Show vertices" is enabled in the View panel

## Example Session

1. Load woody-lo.obj
2. Toggle Deformation Mode ON
3. Click 5-6 vertices at the bottom → Red spheres appear
4. Shift+Click 2-3 vertices at the top → Green spheres appear
5. Switch to "Move vertex" mode
6. Drag green handles upward and sideways
7. Click "Apply Deformation"
8. Watch mesh smoothly deform!
9. (Optional) Add more handles and apply again
10. Use "Reset Mesh" when done

## Advanced Features

### Multiple Iterations

You can apply deformation multiple times:
1. Apply first deformation
2. Keep or modify handles
3. Drag movable handles to new positions
4. Apply again

Each application uses the current mesh state as the starting point.

### Weight Comparison

Try both weight types to see the difference:
- **Cotangent**: Respects triangle geometry, better for irregular meshes
- **Uniform**: Treats all edges equally, simpler but may distort

### Handle Management

- Handles persist across multiple applications
- Use "Clear Handles" to remove all and start fresh
- Individual handles can be toggled on/off by clicking again

## Performance Notes

- **Small meshes** (< 1000 vertices): Near-instant deformation
- **Medium meshes** (1000-10000 vertices): 1-2 seconds
- **Large meshes** (> 10000 vertices): May take several seconds

The solve time is dominated by the sparse LU decomposition, which is O(n^1.5) for typical meshes.

## Limitations

1. **No rotation fitting**: Unlike ARAP, this is purely linear Laplacian deformation
2. **Volume not preserved**: Mesh may shrink or expand
3. **Large rotations**: May appear distorted (use ARAP for large rotations)
4. **Triangle meshes only**: Polygonal faces are skipped

## Output

The deformed mesh is:
- Displayed immediately in the viewer
- Stored in the application (use "Reset Mesh" to undo)
- Can be exported using the provided mesh I/O (not in GUI, use CLI)

---

**Enjoy deforming your meshes!**
