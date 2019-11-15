import math
# Blender
import bpy
import bpy_extras
import mathutils

# define function
def pbrt_disk(height = 0.0, radius = 1.0, innerradius = 0.0, phimax =
              360.0, steps = 32):
    # check arguments (and correct them)
    new_height = height
    if radius <= 0.0:
        print("WARNING: radius <= 0.0 (%s <= 0.0)" % radius)
        new_radius = 1.0
    else:
        new_radius = radius
    if innerradius < 0.0:
        print("WARNING: innerradius < 0.0 (%s < 0.0)" % innerradius)
        new_innerradius = 0.0
    elif innerradius >= new_radius:
        print("WARNING: innerradius >= radius (%s >= %s)" % (innerradius, new_radius))
        new_innerradius = 0.0
    else:
        new_innerradius = innerradius
    if phimax <= 0.0:
        print("WARNING: phimax <= 0.0 (%s <= 0.0)" % phimax)
        new_phimax = 360.0
    elif phimax > 360.0:
        print("WARNING: phimax > 360.0 (%s > 360.0)" % phimax)
        new_phimax = 360.0
    else:
        new_phimax = phimax
    print("pbrt_disk(height = %s, radius = %s, innerradius = %s, phimax = %s)" %
          (new_height, new_radius, new_innerradius, new_phimax))
    # some edges to spin
    verts = [mathutils.Vector((new_innerradius, 0, new_height)),
             mathutils.Vector((new_radius, 0, new_height))]
    edges = [[0, 1]]
    faces = []
    mesh = bpy.data.meshes.new(name="PbrtDisk")
    mesh.from_pydata(verts, edges, faces)
    # mesh.validate(verbose=True)
    bpy_extras.object_utils.object_data_add(bpy.context, mesh)
    # spin
    bpy.ops.object.editmode_toggle()
    bpy.ops.mesh.spin(steps=steps, angle=-math.radians(new_phimax), axis=(0,0,1))
    bpy.ops.mesh.select_all(action='TOGGLE')
    bpy.ops.mesh.select_all(action='TOGGLE')
    if new_phimax == 360.0 or new_innerradius == 0.0:
        bpy.ops.mesh.remove_doubles()
    bpy.ops.object.editmode_toggle()
    # shading: smooth
    bpy.ops.object.shade_smooth()
    # custom properties
    obj = bpy.data.objects["PbrtDisk"]
    obj['_RNA_UI'] = {}
    obj["height"] = new_height
    obj["radius"] = new_radius
    obj["innerradius"] = new_innerradius
    obj["phimax"] = new_phimax

# call function
height = 0.0
radius = 1.0
innerradius = 0.0
phimax = 360.0
steps = 32
pbrt_disk(height = height, radius = radius, innerradius = innerradius,
          phimax = phimax, steps = steps)
