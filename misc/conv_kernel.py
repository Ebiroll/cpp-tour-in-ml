import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np

fig, axes = plt.subplots(ncols=3, figsize=(10, 4))

###############################################################################
# 1. Cartoon Face
###############################################################################
ax_face = axes[0]
ax_face.set_title("Input Image")

# Draw a head (circle)
face_circle = Circle((0.5, 0.5), 0.4, edgecolor='black', facecolor='#ffe0cc')
ax_face.add_patch(face_circle)

# Draw simple hair rectangle on top
ax_face.add_patch(Rectangle((0.15, 0.72), 0.7, 0.15, facecolor='brown', edgecolor='none'))

# Draw eyes (two smaller filled circles)
ax_face.add_patch(Circle((0.35, 0.55), 0.03, color='black'))
ax_face.add_patch(Circle((0.65, 0.55), 0.03, color='black'))

# Draw eyebrows (rectangles)
ax_face.add_patch(Rectangle((0.30, 0.59), 0.10, 0.005, facecolor='black', angle=2))
ax_face.add_patch(Rectangle((0.60, 0.59), 0.10, 0.005, facecolor='black', angle=-2))

# Draw nose (small vertical line)
ax_face.plot([0.50, 0.50], [0.45, 0.50], color='black', linewidth=1.5)

# Draw mouth (simple line)
ax_face.plot([0.40, 0.60], [0.40, 0.40], color='black', linewidth=2)

ax_face.set_xlim([0,1])
ax_face.set_ylim([0,1])
ax_face.set_aspect('equal')
ax_face.axis('off')

###############################################################################
# 2. Convolution Kernel Illustration
###############################################################################
ax_kernel = axes[1]
ax_kernel.set_title("Convolution Kernel")

# We'll just represent the kernel as a small square containing an “eye” shape.
kernel_size = 0.8
kernel_rect = Rectangle((0.1, 0.1), kernel_size, kernel_size,
                        edgecolor='black', facecolor='#ffe0cc')
ax_kernel.add_patch(kernel_rect)

# Eye inside the kernel
ax_kernel.add_patch(Rectangle((0.3, 0.45), 0.4, 0.1, facecolor='black'))  # eyebrow-ish
ax_kernel.add_patch(Circle((0.5, 0.4), 0.08, facecolor='black'))

ax_kernel.set_xlim([0,1])
ax_kernel.set_ylim([0,1])
ax_kernel.set_aspect('equal')
ax_kernel.axis('off')

###############################################################################
# 3. Activation Map (showing bright spots where features match)
###############################################################################
ax_activation = axes[2]
ax_activation.set_title("Activation Map")

# We’ll create a black background
activation_map = np.zeros((100, 100))

# “Activate” two bright spots corresponding to eyes in the input
activation_map[40:50, 30:40] = 1.0  # left eye
activation_map[40:50, 60:70] = 1.0  # right eye

ax_activation.imshow(activation_map, cmap='binary')
ax_activation.axis('off')

plt.tight_layout()
plt.show()

