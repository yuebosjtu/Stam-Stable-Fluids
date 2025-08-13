import os
import argparse
import cv2
from natsort import natsorted

def images_to_video(image_folder, output_path, fps=15, ext='ppm'):
    # 获取所有图片文件
    images = [img for img in os.listdir(image_folder) if img.endswith(f'.{ext}')]
    if not images:
        print(f"No .{ext} images found in {image_folder}")
        return
    images = natsorted(images)
    first_image_path = os.path.join(image_folder, images[0])
    frame = cv2.imread(first_image_path)
    if frame is None:
        print(f"Failed to read {first_image_path}")
        return
    height, width, layers = frame.shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(output_path, fourcc, fps, (width, height))
    for image in images:
        img_path = os.path.join(image_folder, image)
        frame = cv2.imread(img_path)
        if frame is None:
            print(f"Warning: failed to read {img_path}, skipping.")
            continue
        video.write(frame)
    video.release()
    print(f"Video saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description='将文件夹中的图片序列合成为mp4视频')
    parser.add_argument('--input', '-i', required=True, help='图片文件夹路径')
    parser.add_argument('--output', '-o', required=True, help='输出视频文件路径 (如 output.mp4)')
    parser.add_argument('--fps', type=int, default=15, help='帧率 (默认15)')
    parser.add_argument('--ext', default='ppm', help='图片扩展名 (默认ppm)')
    args = parser.parse_args()
    images_to_video(args.input, args.output, args.fps, args.ext)

if __name__ == "__main__":
    main()
