#include <opencv2/opencv.hpp> 
#include <iostream> 
#include <chrono> 
#include <sstream> 

int main() {
    cv::VideoCapture cap(0); 

    if (!cap.isOpened()) { // Проверяем, успешно ли открыт видеопоток
        std::cerr << "Error opening video stream or file" << std::endl;
        return -1; 
    }

    auto start_total = std::chrono::high_resolution_clock::now(); // Запоминаем начальное время для расчета общего времени работы
    auto start = std::chrono::high_resolution_clock::now(); // Время работы программы
    int frame_count = 0; // Счетчик обработанных кадров
    long long total_input_time = 0; // Общее время на чтение кадров
    long long total_transform_time = 0; // Общее время на обработку кадров
    long long total_display_time = 0;  // Общее время на отображение кадров


    while (true) { 
        auto start_input = std::chrono::high_resolution_clock::now(); // Запоминаем время начала чтения кадра
        cv::Mat frame; // Создаем переменную для хранения кадра
        bool ret = cap.read(frame); // Читаем кадр из видеопотока
        auto end_input = std::chrono::high_resolution_clock::now(); // Запоминаем время окончания чтения кадра
        total_input_time += std::chrono::duration_cast<std::chrono::microseconds>(end_input - start_input).count(); 

        if (!ret) { // Проверяем, успешно ли считался кадр
            break;
        }

        auto start_transform = std::chrono::high_resolution_clock::now(); // Запоминаем время начала обработки кадра
        cv::bitwise_not(frame, frame); // Инвертируем цвета кадра

        std::stringstream ss; 
        auto end_time = std::chrono::high_resolution_clock::now(); // Запоминаем текущее время для расчета FPS
        auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_total); // Вычисляем время, прошедшее с начала работы
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start);

        double fps = (double)frame_count / (elapsed_time.count() / 1000.0); 
        ss << "FPS: " << fps; 

        cv::putText(frame, ss.str(), cv::Point(10, 30), cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(0, 255, 0), 2); 
        auto end_transform = std::chrono::high_resolution_clock::now(); 
        total_transform_time += std::chrono::duration_cast<std::chrono::microseconds>(end_transform - start_transform).count(); 
     

        if (elapsed_time.count() >= 1000) { // Выводим статистику каждые 1000 мс (1 секунда)
            long long total_time = total_input_time + total_transform_time + total_display_time; // Общее время обработки

            std::cout << "All time: " << (double)duration.count()/ 1000.0 << " s" << std::endl;
            std::cout << "Input Time: " << (double)total_input_time / 1000.0 << "ms(" << (double)total_input_time * 100.0 / total_time << "%)" << std::endl;
            std::cout << "Transform Time: " << (double)total_transform_time / 1000.0 << "ms(" << (double)total_transform_time * 100.0 / total_time << "%)" << std::endl;
            std::cout << "Display Time: " << (double)total_display_time / 1000.0 << "ms(" << (double)total_display_time * 100.0 / total_time << "%)" << std::endl;
            std::cout << "Total Time: " << (double)total_time / 1000.0 << "ms\n" << std::endl;

            start_total = std::chrono::high_resolution_clock::now(); // Сбрасываем начальное время для следующей секунды
            frame_count = 0; 
            total_input_time = 0; 
            total_transform_time = 0; 
            total_display_time = 0; 
        }

        auto start_display = std::chrono::high_resolution_clock::now(); // Время работы дисплея
        cv::imshow("Video", frame); // Отображаем кадр в окне "Video"
        auto end_display = std::chrono::high_resolution_clock::now();
        total_display_time += std::chrono::duration_cast<std::chrono::microseconds>(end_display - start_display).count(); 

        frame_count++; 
        if (cv::waitKey(1) == ' ') { 
            break; 
        }
    }

    cap.release(); // Закрываем видеопоток
    cv::destroyAllWindows(); // Закрываем все окна OpenCV
    return 0; 
}
