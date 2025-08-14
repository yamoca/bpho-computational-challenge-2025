import sys
import math
import pygame
from pygame import Surface
from PIL import Image
import numpy as np

# ----------------------------- config -----------------------------
WIDTH, HEIGHT = 1000, 700
BG = (245, 246, 250)
AXIS = (30, 30, 30)
LENS = (40, 100, 255)
FOCAL = (40, 100, 255)
OBJECT_BOX = (0, 0, 0)
TEXT = (10, 10, 10)
REAL_CLR = (0, 150, 0)
VIRT_CLR = (230, 120, 20)

# world-to-screen scaling (pixels per "meter" unit)
SCALE = 220
# screen origin at center
ORIGIN = (WIDTH // 2, HEIGHT // 2)

HELP_TEXT = [
    "Controls:",
    "  A/D: object x  |  W/S: object y",
    "  Q/E: focal length f",
    "  Z/X: object width",
    "  R: reset",
    "  Drag & drop an image file into the window",
]

# ------------------------- helpers -------------------------
def world_to_screen(wx, wy):
    cx, cy = ORIGIN
    return (int(cx + wx * SCALE), int(cy - wy * SCALE))

def draw_axis(surface):
    # x-axis
    x0, y0 = world_to_screen(-WIDTH, 0)
    x1, y1 = world_to_screen(WIDTH, 0)
    pygame.draw.line(surface, AXIS, (x0, y0), (x1, y1), 1)

def draw_lens(surface):
    # vertical line at x=0
    (sx, _top) = world_to_screen(0, 999)
    (_, sy) = world_to_screen(0, -999)
    pygame.draw.line(surface, LENS, (sx, _top), (sx, sy), 4)

def draw_focal_points(surface, f, font):
    # +f and -f on x-axis
    for fx in (f, -f):
        px, py = world_to_screen(fx, 0)
        pygame.draw.circle(surface, FOCAL, (px, py), 5)
        label = font.render("F", True, FOCAL)
        surface.blit(label, (px - label.get_width() // 2, py + 8))

def pil_to_surface(pil_img):
    """Convert PIL image to Pygame Surface."""
    mode = pil_img.mode
    size = pil_img.size
    data = pil_img.tobytes()
    return pygame.image.fromstring(data, size, mode).convert_alpha()

def load_image(path):
    try:
        im = Image.open(path).convert("RGBA")
        return im
    except Exception as e:
        print(f"Failed to load image: {e}")
        return None

def blit_image_world(surface, img_surface, left, right, bottom, top):
    """Blit an image to a rectangle given by world coords [left,right]x[bottom,top]."""
    # Compute screen rect
    sx_left, sy_top = world_to_screen(left, top)
    sx_right, sy_bottom = world_to_screen(right, bottom)
    w = max(1, sx_right - sx_left)
    h = max(1, sy_bottom - sy_top)
    scaled = pygame.transform.smoothscale(img_surface, (w, h))
    surface.blit(scaled, (sx_left, sy_top))

def draw_text_block(surface, lines, x, y, font, color=TEXT):
    for i, line in enumerate(lines):
        s = font.render(line, True, color)
        surface.blit(s, (x, y + i * (s.get_height() + 2)))

# ----------------- thin-lens mapping (your logic) -----------------
def lens_image_position(x, y, f):
    """
    Given object at (x, y) with lens at x=0 and focal length f (>0),
    return (is_real, X, Y) or (None, None, None) if x == f (image at infinity).
    """
    if abs(x - f) < 1e-9:
        return None, None, None  # no finite image
    if x > f:
        # real image (inverted)
        X = -f * x / (x - f)
        Y = (y / x) * X
        return True, X, Y
    else:
        # virtual image (upright)
        X = f * x / (x - f)
        Y = (y / x) * X
        return False, X, Y

# ----------------------------- main -----------------------------
def main():
    pygame.init()
    pygame.display.set_caption("Thin Lens – Pygame (simple)")
    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    clock = pygame.time.Clock()
    font = pygame.font.SysFont(None, 20)
    bigfont = pygame.font.SysFont(None, 24)

    # defaults
    f = 1.0
    obj_x = 0.8
    obj_y = 0.0
    obj_w = 0.5

    # load image (arg or placeholder)
    pil_img = None
    if len(sys.argv) > 1:
        pil_img = load_image(sys.argv[1])

    if pil_img is None:
        # simple placeholder image (colored rect with diagonal)
        placeholder = Surface((200, 200), pygame.SRCALPHA)
        placeholder.fill((255, 255, 255, 255))
        pygame.draw.rect(placeholder, OBJECT_BOX, placeholder.get_rect(), 3)
        pygame.draw.line(placeholder, OBJECT_BOX, (0, 0), (200, 200), 2)
        pygame.draw.line(placeholder, OBJECT_BOX, (200, 0), (0, 200), 2)
        img_surface = placeholder
        aspect = 1.0
        img_loaded = False
    else:
        img_surface = pil_to_surface(pil_img)
        aspect = pil_img.height / pil_img.width
        img_loaded = True

    running = True
    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            # drag & drop support
            if event.type == pygame.DROPFILE:
                path = event.file
                newimg = load_image(path)
                if newimg:
                    pil = newimg
                    img_surface = pil_to_surface(pil)
                    aspect = pil.height / pil.width
                    img_loaded = True

            # keys
            if event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False
                # position
                elif event.key == pygame.K_a:  # left
                    obj_x = max(0.05, obj_x - 0.05)
                elif event.key == pygame.K_d:  # right
                    obj_x = min(3.0, obj_x + 0.05)
                elif event.key == pygame.K_w:  # up
                    obj_y = min(1.5, obj_y + 0.05)
                elif event.key == pygame.K_s:  # down
                    obj_y = max(-1.5, obj_y - 0.05)
                # focal length
                elif event.key == pygame.K_q:
                    f = max(0.1, round(f - 0.05, 3))
                elif event.key == pygame.K_e:
                    f = min(2.0, round(f + 0.05, 3))
                # width
                elif event.key == pygame.K_z:
                    obj_w = max(0.2, round(obj_w - 0.05, 3))
                elif event.key == pygame.K_x:
                    obj_w = min(1.0, round(obj_w + 0.05, 3))
                # reset
                elif event.key == pygame.K_r:
                    f, obj_x, obj_y, obj_w = 1.0, 0.8, 0.0, 0.5

        # compute object and image rectangles (world coords)
        obj_h = aspect * obj_w
        obj_left, obj_right = obj_x, obj_x + obj_w
        obj_bottom, obj_top = obj_y, obj_y + obj_h

        # image mapping
        is_real, X, Y = lens_image_position(obj_x, obj_y, f)

        screen.fill(BG)
        draw_axis(screen)
        draw_lens(screen)
        draw_focal_points(screen, f, font)

        # draw object
        blit_image_world(screen, img_surface, obj_left, obj_right, obj_bottom, obj_top)

        # draw image if defined
        status_lines = [
            f"f = {f:.2f}   x = {obj_x:.2f}   y = {obj_y:.2f}   width = {obj_w:.2f}"
        ]

        if is_real is None:
            status_lines.append("Image: none (object at focal point)")
        else:
            if is_real:
                # real (inverted): flip vertically
                img_draw = pygame.transform.flip(img_surface, False, True)
                img_left, img_right = X - obj_w, X
                img_bottom, img_top = Y - obj_h, Y
                blit_image_world(screen, img_draw, img_left, img_right, img_bottom, img_top)
                tag = bigfont.render("Real Image", True, REAL_CLR)
            else:
                # virtual (upright)
                img_left, img_right = X, X + obj_w
                img_bottom, img_top = Y, Y + obj_h
                blit_image_world(screen, img_surface, img_left, img_right, img_bottom, img_top)
                tag = bigfont.render("Virtual Image", True, VIRT_CLR)

            # tag near image
            if is_real is not None:
                tx, ty = world_to_screen(X, Y)
                screen.blit(tag, (tx - tag.get_width() // 2, ty - 28))

        # labels and help
        draw_text_block(screen, status_lines, 12, 12, font)
        draw_text_block(screen, HELP_TEXT, 12, HEIGHT - 16 * (len(HELP_TEXT) + 1), font)

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()

if __name__ == "__main__":
    main()