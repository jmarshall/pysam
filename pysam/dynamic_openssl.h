#ifndef DYNAMIC_OPENSSL_H
#define DYNAMIC_OPENSSL_H

#include <openssl/hmac.h>
#include <openssl/sha.h>

#ifdef DYNAMIC_NETWORK_LIBS

struct dynamic_libcrypto {
    int opened;
    void *libcrypto;
    const EVP_MD *(*f_EVP_sha1)();
    const EVP_MD *(*f_EVP_sha256)();
    unsigned char *(*f_HMAC)(const EVP_MD *, const void *, int, const unsigned char *, size_t, unsigned char *, unsigned int *);
    unsigned char *(*f_SHA256)(const unsigned char *, size_t, unsigned char *);
};

extern struct dynamic_libcrypto libcrypto_table;

#define EVP_sha1   libcrypto_table.f_EVP_sha1
#define EVP_sha256 libcrypto_table.f_EVP_sha256
#define HMAC       libcrypto_table.f_HMAC
#define SHA256     libcrypto_table.f_SHA256

int dynamic_crypto_init();
void dynamic_crypto_cleanup();

#else

#define dynamic_crypto_init()  0
#define dynamic_crypto_cleanup()

#endif
#endif
