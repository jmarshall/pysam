#ifndef DYNAMIC_CURL_H
#define DYNAMIC_CURL_H

#include <curl/curl.h>

#ifdef DYNAMIC_NETWORK_LIBS

struct dynamic_libcurl {
    int opened;
    void *libcurl;
    void (*f_curl_easy_cleanup)(CURL *);
    CURL *(*f_curl_easy_duphandle)(CURL *);
    CURLcode (*f_curl_easy_getinfo)(CURL *, CURLINFO, ...);
    CURL *(*f_curl_easy_init)();
    CURLcode (*f_curl_easy_pause)(CURL *, int);
    CURLcode (*f_curl_easy_perform)(CURL *);
    void (*f_curl_easy_reset)(CURL *);
    CURLcode (*f_curl_easy_setopt)(CURL *, CURLoption, ...);
    const char *(*f_curl_easy_strerror)(CURLcode);
    CURLMcode (*f_curl_multi_add_handle)(CURLM *, CURL *);
    CURLMcode (*f_curl_multi_cleanup)(CURLM *);
    CURLMcode (*f_curl_multi_fdset)(CURLM *, fd_set *, fd_set *, fd_set *, int *);
    CURLMsg *(*f_curl_multi_info_read)(CURLM *, int *);
    CURLM *(*f_curl_multi_init)();
    CURLMcode (*f_curl_multi_perform)(CURLM *, int *);
    CURLMcode (*f_curl_multi_remove_handle)(CURLM *, CURL *);
    const char *(*f_curl_multi_strerror)(CURLMcode);
    CURLMcode (*f_curl_multi_timeout)(CURLM *, long *);
    CURLSHcode (*f_curl_share_cleanup)(CURLSH *);
    CURLSH *(*f_curl_share_init)();
    CURLSHcode (*f_curl_share_setopt)(CURLSH *, CURLSHoption, ...);
    struct curl_slist *(*f_curl_slist_append)(struct curl_slist *, const char *);
    void (*f_curl_slist_free_all)(struct curl_slist *);
    curl_version_info_data *(*f_curl_version_info)(CURLversion);
};

extern struct dynamic_libcurl libcurl_table;

#define curl_easy_cleanup        libcurl_table.f_curl_easy_cleanup
#define curl_easy_duphandle      libcurl_table.f_curl_easy_duphandle
#undef  curl_easy_getinfo
#define curl_easy_getinfo        libcurl_table.f_curl_easy_getinfo
#define curl_easy_init           libcurl_table.f_curl_easy_init
#define curl_easy_pause          libcurl_table.f_curl_easy_pause
#define curl_easy_perform        libcurl_table.f_curl_easy_perform
#define curl_easy_reset          libcurl_table.f_curl_easy_reset
#undef  curl_easy_setopt
#define curl_easy_setopt         libcurl_table.f_curl_easy_setopt
#define curl_easy_strerror       libcurl_table.f_curl_easy_strerror
#define curl_multi_add_handle    libcurl_table.f_curl_multi_add_handle
#define curl_multi_cleanup       libcurl_table.f_curl_multi_cleanup
#define curl_multi_fdset         libcurl_table.f_curl_multi_fdset
#define curl_multi_info_read     libcurl_table.f_curl_multi_info_read
#define curl_multi_init          libcurl_table.f_curl_multi_init
#define curl_multi_perform       libcurl_table.f_curl_multi_perform
#define curl_multi_remove_handle libcurl_table.f_curl_multi_remove_handle
#define curl_multi_strerror      libcurl_table.f_curl_multi_strerror
#define curl_multi_timeout       libcurl_table.f_curl_multi_timeout
#define curl_share_cleanup       libcurl_table.f_curl_share_cleanup
#define curl_share_init          libcurl_table.f_curl_share_init
#undef  curl_share_setopt
#define curl_share_setopt        libcurl_table.f_curl_share_setopt
#define curl_slist_append        libcurl_table.f_curl_slist_append
#define curl_slist_free_all      libcurl_table.f_curl_slist_free_all
#define curl_version_info        libcurl_table.f_curl_version_info

CURLcode dynamic_curl_global_init(long flags);
void dynamic_curl_global_cleanup();

#define curl_global_init         dynamic_curl_global_init
#define curl_global_cleanup      dynamic_curl_global_cleanup

#endif
#endif
