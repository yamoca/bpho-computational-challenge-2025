import pygame
import sys

pygame.init()

SCREEN_WIDTH = 800
SCREEN_HEIGHT = 600
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
pygame.display.set_caption("Interactive Plane Mirror Simulation (Pixel Transform)")

WHITE = (255, 255, 255)
BLUE = (0, 0, 255)

# transformation function
def create_reflected_image(original_image):
    #Creates a new Pygame Surface that is a horizontally flipped version
    #of the original image by transforming each pixel's coordinate.
    width, height = original_image.get_size()

    reflected_image = pygame.Surface((width, height), pygame.SRCALPHA)

    original_image.lock() # lock for faster operations
    reflected_image.lock()

    for y in range(height):
        for x in range(width):
            color = original_image.get_at((x, y))
            # reflect x coordinate
            new_x = width - 1 - x
            reflected_image.set_at((new_x, y), color)

    original_image.unlock()
    reflected_image.unlock()
    
    return reflected_image


try:
    original_object_image = pygame.image.load('object_image.png').convert_alpha()
    object_image = pygame.transform.scale(original_object_image, (120, 150))
except pygame.error as e:
    print(f"Error loading image: {e}")
    print("Make sure 'object_image.png' is in the same directory.")
    pygame.quit()
    sys.exit()


virtual_image = create_reflected_image(object_image)

#init position
object_rect = object_image.get_rect(center=(SCREEN_WIDTH * 0.75, SCREEN_HEIGHT / 2))

# position mirror in middle of screen
mirror_x = SCREEN_WIDTH // 2

# main loop
running = True
dragging = False

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        
        elif event.type == pygame.MOUSEBUTTONDOWN:
            if object_rect.collidepoint(event.pos):
                dragging = True
                mouse_x, mouse_y = event.pos
                offset_x = object_rect.x - mouse_x
                offset_y = object_rect.y - mouse_y
        
        elif event.type == pygame.MOUSEBUTTONUP:
            dragging = False
            
        elif event.type == pygame.MOUSEMOTION:
            if dragging:
                mouse_x, mouse_y = event.pos
                object_rect.x = mouse_x + offset_x
                object_rect.y = mouse_y + offset_y
    
    
    # Calculate the position of the virtual image
    distance_to_mirror = object_rect.centerx - mirror_x
    virtual_image_x = mirror_x - distance_to_mirror
    virtual_image_y = object_rect.centery
    
    virtual_image_rect = virtual_image.get_rect(center=(virtual_image_x, virtual_image_y))

    #drawing to screen  
    # clear the screen
    screen.fill(WHITE)

    # draw mirror 
    pygame.draw.line(screen, BLUE, (mirror_x, 0), (mirror_x, SCREEN_HEIGHT), 3)
    
    screen.blit(object_image, object_rect)
    screen.blit(virtual_image, virtual_image_rect)

    # update the display
    # side note flip seems like a funny name for an update function but its because of double buffering 
    pygame.display.flip()


pygame.quit()
sys.exit()