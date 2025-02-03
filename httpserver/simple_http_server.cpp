// examples / http - server.cpp - *-C++ - *-
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#include <beman/net/net.hpp>
#include <beman/execution/execution.hpp>
#include "demo_algorithm.hpp"
#include "demo_scope.hpp"
#include "demo_task.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string_view>
#include <unordered_map>

namespace ex = beman::execution;
namespace net = beman::net;
using namespace std::chrono_literals;

// ----------------------------------------------------------------------------

std::unordered_map<std::string, std::string> files{
    {"/", "_deps/net-src/examples/data/index.html"},
    {"/favicon.ico", "_deps/net-src/examples/data/favicon.ico"},
    {"/logo.png", "_deps/net-src/examples/data/logo.png"},
};

auto process_request(auto &stream, std::string request) -> demo::task<>
{
    std::istringstream in(request);
    std::string method, url, version;
    if (!(in >> method >> url >> version) || method != "GET")
    {
        std::cout << "not a [supported] HTTP request\n";
        co_return;
    }
    auto it(files.find(url));
    std::cout << "url='" << url << "' -> " << (it == files.end() ? "not found" : it->second) << "\n";
    std::ostringstream out;
    out << std::ifstream(it == files.end() ? std::string() : it->second).rdbuf();
    auto body{out.str()};
    out.clear();
    out.str(std::string());
    out << "HTTP/1.1 " << (it == files.end() ? "404 not found" : "200 found\r\n")
        << "Content-Length: " << body.size() << "\r\n\r\n"
        << body;
    auto response(out.str());
    co_await net::async_send(stream, net::buffer(response));
}

auto timeout(auto scheduler, auto duration, auto sender)
{
    return demo::when_any(
        std::move(sender),
        net::resume_after(scheduler, duration) | demo::into_error([]
                                                                  { return std::error_code(); }));
}

auto make_client(auto scheduler, auto stream) -> demo::task<>
{
    char buffer[16];
    std::string request;
    try
    {
        while (auto n = co_await timeout(scheduler, 3s, net::async_receive(stream, net::buffer(buffer))))
        {
            std::string_view sv(buffer, n);
            request += sv;
            if (request.npos != sv.find("\r\n\r\n"))
            {
                co_await process_request(stream, std::move(request));
                break;
            }
        }
    }
    catch (...)
    {
        std::cout << "ex: timeout\n";
    }
    std::cout << "client done\n";
}

auto main() -> int
{
    net::io_context context;
    demo::scope scope;
    net::ip::tcp::endpoint ep(net::ip::address_v4::any(), 12345);
    net::ip::tcp::acceptor server(context, ep);
    std::cout << "listening on " << ep << "\n";

    scope.spawn(std::invoke(
        [](auto scheduler, auto &scp, auto &svr) -> demo::task<>
        {
            while (true)
            {
                auto [stream, address] = co_await net::async_accept(svr);
                std::cout << "received connection from " << address << "\n";
                scp.spawn(make_client(scheduler, std::move(stream)));
            }
        },
        context.get_scheduler(),
        scope,
        server));

    context.run();
}

#if 0
// simple_http_server.cpp
//
// A sample Sender/Receiver-based HTTP server using the Beman Project libraries.
// This code uses Beman::execution (formerly execution26) and Beman::net (formerly net29)
// to asynchronously accept TCP connections and reply with a fixed HTTP response.
//
// Note: This is a simplified example meant for illustration purposes.
// Actual error handling, request parsing, and API details may vary depending on
// the current state of the Beman libraries.
//
// Required headers (provided by Beman Project):
// #include <beman/execution/execution.hpp> // For default_scheduler, just(), connect(), etc. :contentReference[oaicite:2]{index=2}
// #include <beman/net/tcp_listener.hpp>    // For tcp_listener and async_accept() :contentReference[oaicite:3]{index=3}
// #include <beman/net/tcp_socket.hpp>      // For tcp_socket and async_write()

#include <beman/net/net.hpp>
#include <beman/execution/execution.hpp>
#include <iostream>
#include <string>

using namespace beman::execution;
using namespace beman::net;

// Receiver to handle a single accepted connection.
// When a tcp_socket is received, this receiver sends a fixed HTTP response.
struct connection_receiver
{
    // Called when a connection (tcp_socket) is available.
    void set_value(tcp_socket sock)
    {
        // Prepare a simple HTTP response.
        std::string response =
            "HTTP/1.1 200 OK\r\n"
            "Content-Length: 13\r\n"
            "Content-Type: text/plain\r\n"
            "\r\n"
            "Hello, world!";

        // Issue an asynchronous write on the socket.
        // (We assume async_write returns a sender that completes when the write is done.)
        auto write_sender = sock.async_write(response);
        // When writing is complete, close the connection.
        auto close_sender = write_sender.then([sock = std::move(sock)]() mutable
                                              { sock.close(); });
        // Connect the close operation to a no-op receiver.
        connect(close_sender, []() {});
    }
    void set_error(std::exception_ptr eptr)
    {
        try
        {
            if (eptr)
                std::rethrow_exception(eptr);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Connection error: " << e.what() << "\n";
        }
    }
    void set_stopped()
    {
        // No special action on stop.
    }
};

// Receiver to handle accepted connections from the listener.
// It recursively reissues an asynchronous accept for each new connection.
struct accept_receiver
{
    tcp_listener &listener; // Reference to our listening socket.

    // Called when an accepted tcp_socket is available.
    void set_value(tcp_socket sock)
    {
        // For each accepted connection, create a sender that yields the socket
        // and connect it to our connection_receiver.
        auto conn_sender = just(std::move(sock));
        connection_receiver conn_recv;
        connect(conn_sender, conn_recv);

        // Reissue an asynchronous accept to handle the next incoming connection.
        auto next_accept = listener.async_accept();
        connect(next_accept, *this);
    }
    void set_error(std::exception_ptr eptr)
    {
        try
        {
            if (eptr)
                std::rethrow_exception(eptr);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Accept error: " << e.what() << "\n";
        }
    }
    void set_stopped()
    {
        // No special action on stop.
    }
};

int main()
{
    try
    {
        // Create a default scheduler for asynchronous operations.
        auto sched = default_scheduler();

        // Create a TCP listener on all interfaces ("0.0.0.0") at port 8080.
        // Pass in the scheduler so that the listener uses it for asynchronous I/O.
        tcp_listener listener("0.0.0.0", 8080, sched);

        // Start accepting connections asynchronously.
        // async_accept() returns a sender that will eventually produce a tcp_socket.
        auto accept_sender = listener.async_accept();

        // Create an accept_receiver and connect it to the accept_sender.
        accept_receiver acc_recv{listener};
        connect(accept_sender, acc_recv);

        // Run the scheduler's event loop.
        // This call blocks until no more work remains.
        sched.run();
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error in main: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
#endif