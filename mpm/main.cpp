#include <iostream>
#include <windows.h>
#include <time.h>
#include <stdio.h>
#include "simulator.h"
#include "guicon.h"
#include <crtdbg.h>

using namespace std;

static Simulator *simulator = nullptr;

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);
void DrawPixels(HWND hwnd, BITMAPINFO *);

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
    PWSTR lpCmdLine, int nCmdShow) {

    RedirectIOToConsole();

    if (!simulator) simulator = new Simulator(64, 64 * 64, 1.f, 0.5f);

    MSG  msg;
    WNDCLASSW wc = { 0 };

    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpszClassName = L"Pixels";
    wc.hInstance = hInstance;
    wc.hbrBackground = GetSysColorBrush(COLOR_3DFACE);
    wc.lpfnWndProc = WndProc;
    wc.hCursor = LoadCursor(0, IDC_ARROW);

    RegisterClassW(&wc);
    CreateWindowW(wc.lpszClassName, L"Pixels",
        WS_OVERLAPPEDWINDOW | WS_VISIBLE,
        100, 100, 800, 820, NULL, NULL, hInstance, NULL);


    while (GetMessage(&msg, NULL, 0, 0)) {

        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    srand(time(NULL));

    return (int)msg.wParam;
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg,
    WPARAM wParam, LPARAM lParam) {

    BITMAPINFO bmi;
    ZeroMemory(&bmi, sizeof(BITMAPINFO));
    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biPlanes = 1;
    bmi.bmiHeader.biBitCount = 24;
    bmi.bmiHeader.biCompression = BI_RGB;
    bmi.bmiHeader.biWidth = 800;
    bmi.bmiHeader.biHeight = 800;

    switch (msg) {

    case WM_CREATE:
        SetTimer(hwnd, 1, 1, NULL);
        break;

    case WM_PAINT:
        DrawPixels(hwnd, &bmi);
        break;

    case WM_TIMER:
        InvalidateRect(hwnd, NULL, FALSE);
        break;

    case WM_DESTROY:

        PostQuitMessage(0);
        return 0;
    }

    return DefWindowProcW(hwnd, msg, wParam, lParam);
}

void DrawPixels(HWND hwnd, BITMAPINFO *bmi) {

    PAINTSTRUCT ps;
    RECT r;

    GetClientRect(hwnd, &r);

    if (r.bottom == 0) {

        return;
    }

    HDC hdc = BeginPaint(hwnd, &ps);

    simulator->step();
    Mat<unsigned char> image(800 * 3, 800, fill::zeros);
    simulator->render(image);

    /*for (size_t x = 0; x < 800; x++) {
        for (size_t y = 0; y < 800; y++) {
            if (image(y, x) > 0) SetPixel(hdc, x, y, RGB(image(y, x), 0, 0));
        }
    }*/

    SetDIBitsToDevice(hdc, 0, 0, image.n_rows / 3, image.n_cols, 0, 0, 0, image.n_cols, image.memptr(), bmi, DIB_RGB_COLORS);

    /*for (int i = 0; i < 1000; i++) {

        int x = rand() % r.right;
        int y = rand() % r.bottom;
        SetPixel(hdc, x, y, RGB(255, 0, 0));
    }*/

    EndPaint(hwnd, &ps);

}