import math
# Blender
import bpy
import bpy_extras
import mathutils

# define function
def pbrt_cylinder(radius = 1.0, zmin = -1.0, zmax = 1.0, phimax = 360.0, steps = 32):
    # check arguments (and correct them)
    if radius <= 0.0:
        print("WARNING: radius <= 0.0 (%s <= 0.0)" % radius)
        new_radius = 1.0
    else:
        new_radius = radius
    new_zmin = zmin
    new_zmax = zmax
    if phimax <= 0.0:
        print("WARNING: phimax <= 0.0 (%s <= 0.0)" % phimax)
        new_phimax = 360.0
    elif phimax > 360.0:
        print("WARNING: phimax > 360.0 (%s > 360.0)" % phimax)
        new_phimax = 360.0
    else:
        new_phimax = phimax
    print("pbrt_cylinder(radius = %s, zmin = %s, zmax = %s, phimax = %s)" %
          (new_radius, new_zmin, new_zmax, new_phimax))
    # some edges to spin
    verts = [mathutils.Vector((new_radius, 0, new_zmin)),
             mathutils.Vector((new_radius, 0, new_zmax))]
    edges = [[0, 1]]
    faces = []
    mesh = bpy.data.meshes.new(name="PbrtCylinder")
    mesh.from_pydata(verts, edges, faces)
    mesh.validate(verbose=True)
    bpy_extras.object_utils.object_data_add(bpy.context, mesh)
    # spin
    bpy.ops.object.editmode_toggle()
    bpy.ops.mesh.spin(steps=steps, angle=-math.radians(new_phimax), axis=(0,0,1))
    bpy.ops.mesh.select_all(action='TOGGLE')
    bpy.ops.mesh.select_all(action='TOGGLE')
    bpy.ops.mesh.normals_make_consistent(inside=False)
    if new_phimax == 360.0:
        bpy.ops.mesh.remove_doubles()
    bpy.ops.object.editmode_toggle()
    # shading: smooth
    bpy.ops.object.shade_smooth()
    # custom properties
    obj = bpy.data.objects["PbrtCylinder"]
    obj['_RNA_UI'] = {}
    obj["radius"] = new_radius
    obj["zmin"] = new_zmin
    obj["zmax"] = new_zmax
    obj["phimax"] = new_phimax

# call function
radius = 1.0
zmin = -1.0
zmax = 1.0
phimax = 360.0
steps = 32
pbrt_cylinder(radius = radius, zmin = zmin, zmax = zmax, phimax = phimax, steps = steps)
