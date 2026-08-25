#include <dlfcn.h>

#include "htslib/hts_log.h"

#include "dynamic_curl.h"
#include "dynamic_openssl.h"

struct dynamic_libcurl libcurl_table = { 0, NULL };
static const char *libcurl_soname = NULL;

void dynamic_set_libcurl_soname(const char *soname) {
    libcurl_soname = soname;
}

CURLcode dynamic_curl_global_init(long flags)
{
    if (libcurl_table.opened > 0) {
        libcurl_table.opened++;
        return CURLE_OK;
    }

    const char *tried_sonames;
    void *libcurl;

    if (libcurl_soname) {
        tried_sonames = libcurl_soname;
        libcurl = dlopen(libcurl_soname, RTLD_LAZY | RTLD_LOCAL);
    }
    else {
        tried_sonames = "libcurl.so.4";
        libcurl = dlopen("libcurl.so.4", RTLD_LAZY | RTLD_LOCAL);
    }

    if (libcurl == NULL) {
        hts_log_error("Could not load %s: %s", tried_sonames, dlerror());
        return CURLE_FAILED_INIT;
    }

    CURLcode (*f_curl_global_init)(long);
    *(void **) &f_curl_global_init = dlsym(libcurl, "curl_global_init");
    CURLcode ret = f_curl_global_init? f_curl_global_init(flags) : CURLE_FAILED_INIT;
    if (ret != CURLE_OK) { dlclose(libcurl); return ret; }

    libcurl_table.opened++;
    libcurl_table.libcurl = libcurl;

    *(void **) &libcurl_table.f_curl_easy_cleanup        = dlsym(libcurl, "curl_easy_cleanup");
    *(void **) &libcurl_table.f_curl_easy_duphandle      = dlsym(libcurl, "curl_easy_duphandle");
    *(void **) &libcurl_table.f_curl_easy_getinfo        = dlsym(libcurl, "curl_easy_getinfo");
    *(void **) &libcurl_table.f_curl_easy_init           = dlsym(libcurl, "curl_easy_init");
    *(void **) &libcurl_table.f_curl_easy_pause          = dlsym(libcurl, "curl_easy_pause");
    *(void **) &libcurl_table.f_curl_easy_perform        = dlsym(libcurl, "curl_easy_perform");
    *(void **) &libcurl_table.f_curl_easy_reset          = dlsym(libcurl, "curl_easy_reset");
    *(void **) &libcurl_table.f_curl_easy_setopt         = dlsym(libcurl, "curl_easy_setopt");
    *(void **) &libcurl_table.f_curl_easy_strerror       = dlsym(libcurl, "curl_easy_strerror");
    *(void **) &libcurl_table.f_curl_multi_add_handle    = dlsym(libcurl, "curl_multi_add_handle");
    *(void **) &libcurl_table.f_curl_multi_cleanup       = dlsym(libcurl, "curl_multi_cleanup");
    *(void **) &libcurl_table.f_curl_multi_fdset         = dlsym(libcurl, "curl_multi_fdset");
    *(void **) &libcurl_table.f_curl_multi_info_read     = dlsym(libcurl, "curl_multi_info_read");
    *(void **) &libcurl_table.f_curl_multi_init          = dlsym(libcurl, "curl_multi_init");
    *(void **) &libcurl_table.f_curl_multi_perform       = dlsym(libcurl, "curl_multi_perform");
    *(void **) &libcurl_table.f_curl_multi_remove_handle = dlsym(libcurl, "curl_multi_remove_handle");
    *(void **) &libcurl_table.f_curl_multi_strerror      = dlsym(libcurl, "curl_multi_strerror");
    *(void **) &libcurl_table.f_curl_multi_timeout       = dlsym(libcurl, "curl_multi_timeout");
    *(void **) &libcurl_table.f_curl_share_cleanup       = dlsym(libcurl, "curl_share_cleanup");
    *(void **) &libcurl_table.f_curl_share_init          = dlsym(libcurl, "curl_share_init");
    *(void **) &libcurl_table.f_curl_share_setopt        = dlsym(libcurl, "curl_share_setopt");
    *(void **) &libcurl_table.f_curl_slist_append        = dlsym(libcurl, "curl_slist_append");
    *(void **) &libcurl_table.f_curl_slist_free_all      = dlsym(libcurl, "curl_slist_free_all");
    *(void **) &libcurl_table.f_curl_version_info        = dlsym(libcurl, "curl_version_info");

    return CURLE_OK;
}

void dynamic_curl_global_cleanup()
{
    libcurl_table.opened--;
    if (libcurl_table.opened > 0) return;

    void (*f_curl_global_cleanup)();
    *(void **) &f_curl_global_cleanup = dlsym(libcurl_table.libcurl, "curl_global_cleanup");
    if (f_curl_global_cleanup) f_curl_global_cleanup();

    dlclose(libcurl_table.libcurl);
    libcurl_table.libcurl = NULL;
}


struct dynamic_libcrypto libcrypto_table = { 0, NULL };
static const char *libcrypto_soname = NULL;

void dynamic_set_libcrypto_soname(const char *soname) {
    libcrypto_soname = soname;
}

int dynamic_crypto_init()
{
    if (libcrypto_table.opened > 0) {
        libcrypto_table.opened++;
        return 0;
    }

    const char *tried_sonames;
    void *libcrypto;

    if (libcrypto_soname) {
        tried_sonames = libcrypto_soname;
        libcrypto = dlopen(libcrypto_soname, RTLD_LAZY | RTLD_LOCAL);
    }
    else {
        tried_sonames = "libcrypto.so.{4,3,1.1}";
        libcrypto = dlopen("libcrypto.so.4", RTLD_LAZY | RTLD_LOCAL);
        if (libcrypto == NULL) libcrypto = dlopen("libcrypto.so.3", RTLD_LAZY | RTLD_LOCAL);
        if (libcrypto == NULL) libcrypto = dlopen("libcrypto.so.1.1", RTLD_LAZY | RTLD_LOCAL);
    }

    if (libcrypto == NULL) {
        hts_log_error("Could not load %s: %s", tried_sonames, dlerror());
        return -1;
    }

    libcrypto_table.opened++;
    libcrypto_table.libcrypto = libcrypto;

    *(void **) &libcrypto_table.f_EVP_sha1   = dlsym(libcrypto, "EVP_sha1");
    *(void **) &libcrypto_table.f_EVP_sha256 = dlsym(libcrypto, "EVP_sha256");
    *(void **) &libcrypto_table.f_HMAC       = dlsym(libcrypto, "HMAC");
    *(void **) &libcrypto_table.f_SHA256     = dlsym(libcrypto, "SHA256");

    return 0;
}

void dynamic_crypto_cleanup()
{
    libcrypto_table.opened--;
    if (libcrypto_table.opened > 0) return;

    dlclose(libcrypto_table.libcrypto);
    libcrypto_table.libcrypto = NULL;
}
