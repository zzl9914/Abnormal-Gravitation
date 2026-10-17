// freestanding wasm 构建专用：clang 可能把结构体拷贝/清零优化成对
// memcpy / memset / memmove 的调用，而 -nostdlib 下没有 libc，必须自己提供。
// 原生构建不编译本文件（用系统的）。

extern "C" {

void* memcpy(void* dst, const void* src, unsigned long n) {
    unsigned char* d = (unsigned char*)dst;
    const unsigned char* s = (const unsigned char*)src;
    for (unsigned long i = 0; i < n; i++) d[i] = s[i];
    return dst;
}

void* memset(void* dst, int value, unsigned long n) {
    unsigned char* d = (unsigned char*)dst;
    unsigned char v = (unsigned char)value;
    for (unsigned long i = 0; i < n; i++) d[i] = v;
    return dst;
}

void* memmove(void* dst, const void* src, unsigned long n) {
    unsigned char* d = (unsigned char*)dst;
    const unsigned char* s = (const unsigned char*)src;
    if (d == s) return dst;
    if (d < s) {
        for (unsigned long i = 0; i < n; i++) d[i] = s[i];
    } else {
        for (unsigned long i = n; i > 0; i--) d[i - 1] = s[i - 1];
    }
    return dst;
}

}  // extern "C"
