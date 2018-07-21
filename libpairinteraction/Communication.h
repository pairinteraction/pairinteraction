#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <zmq.h>

#include <boost/format.hpp>

namespace zmq {

struct endpoint {
    static std::string name;
};

/** \brief ZeroMQ error
 *
 * Exception to identify ZeroMQ errors.  The constructor has no
 * arguments as ZeroMQ handles errors using a global variable which is
 * evaluated using zmq_errno() inside the exception.
 *
 * \code
 * try
 * {
 *   // code with zeromq
 * }
 * catch (zmq::error &e)
 * {
 *   // handle exception
 * }
 * \endcode
 */
class error : public std::exception {
private:
    int errnum;

public:
    /** \brief Constructor */
    error() : errnum(zmq_errno()) {}

    char const *what() const noexcept override { return zmq_strerror(errnum); }

    /** \brief %what() function
     *
     * \returns error message
     */
    int num() const noexcept { return errnum; }
};

/** \brief Socket class
 *
 * The socket class exposes parts of the ZeroMQ C API in an
 * object-oriented manner.  Only the functionality we require is
 * implemented.  The socket object implicitly behaves like a pointer
 * to void which is what ZeroMQ expects.  Thus a socket object can be
 * used like a regular %socket pointer with the C API.  For example
 * the `zmq_sendmsg` is not implemented but can be used nevertheless.
 * \code
 * if ( zmq_sendmsg(socket, msg, flags) == -1 )
 *   throw zmq::error();
 * \endcode
 * When using the C API directly, errors have to be checked explicitly!
 *
 * For details see http://api.zeromq.org/.
 */
class socket {
    std::unique_ptr<void, int (*)(void *)> m_socket;
    mutable bool use_cout;

public:
    /** \brief Constructor
     *
     * The constructor is not meant to be called by itself but is rather
     * used in the context object to return a new socket.
     *
     * \param[in] new_socket    pre-allocated ZeroMQ socket
     * \throws zmq::error
     */
    explicit socket(void *new_socket) : m_socket(new_socket, zmq_close), use_cout(false) {
        if (new_socket == nullptr) {
            throw error();
        }
    }

    /** \brief Conversion operator
     *
     * A socket object implicitly behaves like a pointer to void which
     * is what ZeroMQ expects, so that a socket object can be used like
     * a regular %socket pointer with the C API.
     */
    operator void *() const noexcept { return m_socket.get(); }

    /** \brief Accept incoming connections on a socket.
     *
     * The bind() function binds the socket to a local endpoint and
     * then accepts incoming connections on that endpoint.
     *
     * \param[in] endpoint   string with transport protocal and address
     * \throws zmq::error
     */
    void bind(char const *endpoint) const {
        if (endpoint == std::string{"use_cout"}) {
            use_cout = true;
            return;
        }
        if (zmq_bind(*this, endpoint) == -1) {
            throw error();
        }
    }

    /** \brief Create outgoing connection from socket
     *
     * The connect() function connects the socket to an endpoint and
     * then accepts incoming connections on that endpoint.
     *
     * \param[in] endpoint   string with transport protocal and address
     * \throws zmq::error
     */
    void connect(char const *endpoint) const {
        if (endpoint == std::string{"use_cout"}) {
            use_cout = true;
            return;
        }
        if (zmq_connect(*this, endpoint) == -1) {
            throw error();
        }
    }

    /** \brief Send a message part on a socket
     *
     * The send() function shall queue a message created from its
     * printf-style arguments.
     *
     * \param[in] fmt    printf-style string for argument formatting
     * \param[in] args   variadic arguments
     * \returns number of bytes written (without zero terminator)
     */
    template <typename... Args>
    int send(char const *fmt, Args &&... args) {
        boost::format formatter(fmt);

        // This is horrible.  In C++17 we have the binary fold expression
        using expander = int[];
        (void)expander{0, (void(formatter % std::forward<Args>(args)), 0)...};

        std::string tmp = formatter.str();
        int len = tmp.size();

        if (use_cout) {
            std::cout << tmp << std::endl;
            return len;
        }

        if (zmq_send(*this, tmp.c_str(), len + 1, 0) == -1) {
            throw error();
        }
        return len;
    }
};

/** \brief Context class
 *
 * The context class exposes parts of the ZeroMQ C API in an
 * object-oriented manner.  Only the functionality we require is
 * implemented.  The context object implicitly behaves like a pointer
 * to void which is what ZeroMQ expects.  Thus a context object can be
 * used like a regular %context pointer with the C API.  For example
 * the `zmq_ctx_set` is not implemented but can be used nevertheless.
 * \code
 * if ( zmq_ctx_set(context, ZMQ_IPV6, 1) == -1 )
 *   throw zmq::error();
 * \endcode
 * When using the C API directly, errors have to be checked explicitly!
 *
 * For details see http://api.zeromq.org/.
 */
class context {
    std::unique_ptr<void, int (*)(void *)> m_context;

public:
    /** \brief Constructor */
    context() noexcept : m_context(zmq_ctx_new(), zmq_ctx_term) {}

    /** \brief Conversion operator
     *
     * A socket object implicitly behaves like a pointer to void which
     * is what ZeroMQ expects, so that a socket object can be used like
     * a regular %socket pointer with the C API.
     */
    operator void *() const noexcept { return m_context.get(); }

    /** Create ZeroMQ socket
     *
     * The socket() function shall create a ZeroMQ socket within the
     * context and return an opaque handle to the newly created socket.
     *
     * \param[in] type   specifies the socket type
     */
    zmq::socket socket(int type) const {
        /** errors of `zmq_socket` are handled in the zmq::socket constructor */
        return zmq::socket(zmq_socket(*this, type));
    }
};

} // namespace zmq
