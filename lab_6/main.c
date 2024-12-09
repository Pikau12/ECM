#include <iostream>
#include <libusb.h>
#include <stdio.h>

void printdev(libusb_device *dev);

int main() {
    libusb_device **devs;// указатель на указатель на устройство   // используется для получения списка устройств
    libusb_context *ctx = nullptr;  // контекст сессии libusb 
    int r; // для возвращаемых значений
    ssize_t cnt; // число найденных USB-устройств
    ssize_t i; // индексная переменная цикла перебора всех устройств
    // инициализировать библиотеку libusb, открыть сессию работы с libusb

    r = libusb_init(&ctx);
    if (r < 0) {
        fprintf(stderr, "Failed to initialize libusb: %d.\n", r);
        return 1;
    }

    libusb_set_debug(ctx, 3);

    cnt = libusb_get_device_list(ctx, &devs);
    if (cnt < 0) {
        fprintf(stderr, "Error getting device list\n");
        return 1;
    }

    printf("Found devices: %d\n\n", (int)cnt);
    for (i = 0; i < cnt; i++) {
        printdev(devs[i]);
    }

    libusb_free_device_list(devs, 1);
    libusb_exit(ctx);
    return 0;
}

void printdev(libusb_device *dev) {
    struct libusb_device_descriptor desc; // Структура для хранения дескриптора устройства.
    libusb_device_handle *handle; // Дескриптор устройства (для работы с ним).
    int r = libusb_get_device_descriptor(dev, &desc); // Получение дескриптора устройства.
    if (r < 0) {
        fprintf(stderr, "Failed to get device descriptor: %s\n", libusb_error_name(r));
        return;
    }

    printf("Device Class: Ox%02x\n", desc.bDeviceClass);
    printf("Vendor ID: Ox%04x\n", desc.idVendor);
    printf("Product ID: Ox%04x\n", desc.idProduct);

    r = libusb_open(dev, &handle);
    if (r != LIBUSB_SUCCESS) {
        fprintf(stderr, "Failed to open device: %s\n", libusb_error_name(r));
        return;
    }


    if (desc.iSerialNumber > 0) {
        unsigned char serial[256];
        r = libusb_get_string_descriptor_ascii(handle, desc.iSerialNumber, serial, sizeof(serial));
        if (r >= 0) {
            printf("Serial Number: %s\n\n", serial);
        } else {
            fprintf(stderr, "Failed to get Serial Number: %s\n", libusb_error_name(r));
        }
    } else {
        printf("No Serial Number\n\n");
    }

    libusb_close(handle);
}
