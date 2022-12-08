import moviepy.editor as mp

clip = mp.VideoFileClip("foo.gif")
clip.write_videofile("output.mp4")
