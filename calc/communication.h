#pragma once

#include <memory>
#include <cstring>
#include <cstdarg>
#include <zmq.h>


namespace zmq {

class error : public std::exception
{
private:
  int errnum;
public:
  error() : errnum(zmq_errno()) {}

  char const * what() const noexcept override
  {
    return zmq_strerror(errnum);
  }

  int num() const noexcept
  {
    return errnum;
  }
};


class socket
{
  std::unique_ptr<void, std::function<void(void*)>> m_socket;
public:
  socket(void * new_socket)
    : m_socket( new_socket, zmq_close )
  {
    if ( !new_socket )
      throw error();
  }
  operator void *() const { return m_socket.get(); }
  void bind(char const *target) const
  {
    if ( zmq_bind(*this, target) != 0 )
      throw error();
  }
  void connect(char const *target) const
  {
    if ( zmq_connect(*this, target) != 0 )
      throw error();
  }
  void send(char const *fmt, ...) const __attribute__((format (printf, 2, 3)));
  void send(char const *fmt, ...)
  {
    va_list args;
    va_start(args, fmt);
    char *tmp = nullptr;
    vasprintf(&tmp, fmt, args);
    zmq_send(*this, tmp, strlen(tmp)+1, 0);
    free(tmp);
  }
};

class context
{
  std::unique_ptr<void, std::function<int(void*)>> m_context;
public:
  context() : m_context( zmq_ctx_new(), zmq_ctx_term ) {}
  operator void *() const { return m_context.get(); }
  zmq::socket socket(int flags) const { return zmq::socket( zmq_socket(*this, flags) ); }
};

}
