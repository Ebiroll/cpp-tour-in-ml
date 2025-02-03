#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <atomic>
#include <websocket.h>
#include <stdexec/execution.hpp>
#include <exec/static_thread_pool.hpp>
#include <cmath> // For std::sqrt

// Include the provided WebSocket code here (omitted for brevity)

struct Body {
    double m;
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
};


std::vector<Body> bodies_;

class Simulation {
public:


Simulation() {
    const double G = 6.67430e-11;  // Gravitational constant
    const double AU = 1.496e11;    // 1 astronomical unit in meters
    const double mass_base = 1e28; // Base mass (similar to small stars)
    
    // Common velocity magnitude calculation
    auto calc_velocity = [&](double m1, double m2, double distance) {
        return std::sqrt(G * (m1 + m2) / distance);
    };

    // Triple star system with chaotic initial conditions
    Body star1, star2, star3;
    
    // Mass configuration (similar but not identical)
    star1.m = mass_base * 1.0;
    star2.m = mass_base * 0.9;
    star3.m = mass_base * 1.1;

    // Initial positions (triangular configuration)
    star1.x =  AU * 0.5;
    star1.y = -AU * 0.3;
    star2.x = -AU * 0.4;
    star2.y =  AU * 0.6;
    star3.x =  AU * 0.2;
    star3.y =  AU * 0.1;

    // Velocity configuration (non-symmetric, non-orbital)
    const double speed_base = 1e3; // 10 km/s base speed
    star1.vx =  speed_base * 0.7;
    star1.vy = -speed_base * 1.1;
    star2.vx = -speed_base * 0.9;
    star2.vy =  speed_base * 0.8;
    star3.vx =  speed_base * 1.2;
    star3.vy =  speed_base * 0.4;

    // Center-of-mass correction
    const double total_mass = star1.m + star2.m + star3.m;
    const double com_vx = (star1.m*star1.vx + star2.m*star2.vx + star3.m*star3.vx) / total_mass;
    const double com_vy = (star1.m*star1.vy + star2.m*star2.vy + star3.m*star3.vy) / total_mass;
    
    // Adjust velocities to maintain system in view
    star1.vx -= com_vx;
    star1.vy -= com_vy;
    star2.vx -= com_vx;
    star2.vy -= com_vy;
    star3.vx -= com_vx;
    star3.vy -= com_vy;

    bodies_.push_back(star1);
    bodies_.push_back(star2);
    bodies_.push_back(star3);
}

# if 0    

Simulation() {
    const double G = 6.67430e-11;  // Gravitational constant
    const double AU = 1.496e11;    // Astronomical unit in meters
    
    // Create Sun
    Body sun;
    sun.m = 1.989e30;
    sun.x = 0;
    sun.y = 0;
    bodies_.push_back(sun);

    // Helper function to create planets
    auto createPlanet = [&](double mass, double au_distance, double orbital_speed) {
        Body planet;
        planet.m = mass;
        planet.x = au_distance * AU;  // Place along x-axis
        planet.vy = orbital_speed;    // Initial tangential velocity
        return planet;
    };

    // Add inner planets (approximate values)
    bodies_.push_back(createPlanet(
        3.301e23,   // Mercury mass
        0.387,      // AU distance
        std::sqrt(G * sun.m / (0.387 * AU))  // Orbital velocity
    ));

    bodies_.push_back(createPlanet(
        4.867e24,   // Venus mass
        0.723,      // AU distance
        std::sqrt(G * sun.m / (0.723 * AU))  // ~35,000 m/s
    ));

    // Add outer planets
    bodies_.push_back(createPlanet(
        5.972e24,   // Earth mass
        1.0,        // AU distance
        std::sqrt(G * sun.m / (1.0 * AU))    // ~29,780 m/s
    ));

    bodies_.push_back(createPlanet(
        6.417e23,   // Mars mass
        1.524,      // AU distance
        std::sqrt(G * sun.m / (1.524 * AU))  // ~24,100 m/s
    ));

    // Add outer gas giants
    bodies_.push_back(createPlanet(
        1.898e27,   // Jupiter mass
        5.203,      // AU distance
        std::sqrt(G * sun.m / (5.203 * AU))  // ~13,060 m/s
    ));

    bodies_.push_back(createPlanet(
        5.683e26,   // Saturn mass
        9.537,      // AU distance
        std::sqrt(G * sun.m / (9.537 * AU))  // ~9,680 m/s
    ));
}
    Simulation() {
        Body sun;
        sun.m = 1.989e30;
        sun.x = 0;
        bodies_.push_back(sun);

        Body earth;
        earth.m = 5.972e24;
        earth.x = 1.496e11; // ~1 AU
        earth.vy = 2.978e4; // Orbital velocity
        bodies_.push_back(earth);
    };
#endif

    void addBody(const Body& body) {
        std::lock_guard<std::mutex> lock(mutex_);
        bodies_.push_back(body);
    }

    void removeBody(int id) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (id >= 0 && id < static_cast<int>(bodies_.size())) {
            bodies_.erase(bodies_.begin() + id);
        }
    }

    void removeMostDistant() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (bodies_.size() > 1) {
            double max_dist = 0;
            size_t max_id = 0;
            for (size_t i = 1; i < bodies_.size(); ++i) {
                const auto& body = bodies_[i];
                double dist = std::sqrt(body.x*body.x + body.y*body.y + body.z*body.z);
                if (dist > max_dist) {
                    max_dist = dist;
                    max_id = i;
                }
            }
            bodies_.erase(bodies_.begin() + max_id);
        }

    }


    void forceBody(int id, double x, double y, double z) {
        #if 0
        std::lock_guard<std::mutex> lock(mutex_);
        if (id >= 0 && id < static_cast<int>(bodies_.size())) {
            auto& body = bodies_[id];
            body.x = x;
            body.y = y;
            body.z = z;
        }
        #endif
    }

    void step(double dt) {
        std::lock_guard<std::mutex> lock(mutex_);
        const double G = 6.67430e-11;
        
        // Reset accelerations
        for (auto& body : bodies_) {
            body.ax = body.ay = body.az = 0.0;
        }

        // Calculate forces
        for (size_t i = 0; i < bodies_.size(); ++i) {
            for (size_t j = i + 1; j < bodies_.size(); ++j) {
                Body& a = bodies_[i];
                Body& b = bodies_[j];
                double dx = b.x - a.x;
                double dy = b.y - a.y;
                double dz = b.z - a.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                // Avoid division by zero with softening
                double F = G * a.m * b.m / (dist*dist + 1e-10);  
                
                // Force components (direction from a to b)
                double Fx = F * dx / dist;
                double Fy = F * dy / dist;
                double Fz = F * dz / dist;

                // Apply equal and opposite forces
                a.ax += Fx / a.m;  // Correct: a is pulled toward b
                a.ay += Fy / a.m;
                a.az += Fz / a.m;
                
                b.ax -= Fx / b.m;  // Correct: b is pulled toward a
                b.ay -= Fy / b.m;
                b.az -= Fz / b.m;
            }
        }

        // Update velocities and positions
        for (auto& body : bodies_) {
            body.vx += body.ax * dt;
            body.vy += body.ay * dt;
            body.vz += body.az * dt;
            body.x += body.vx * dt + 0.5 * body.ax * dt * dt;
            body.y += body.vy * dt + 0.5 * body.ay * dt * dt;
            body.z += body.vz * dt + 0.5 * body.az * dt * dt;
        }
    }
    std::vector<Body> getBodies() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return bodies_;
    }

private:
    mutable std::mutex mutex_;
};


std::string serializeBodies(const std::vector<Body>& bodies) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < bodies.size(); ++i) {
        const auto& body = bodies[i];
        oss << "{"
            << "\"m\":" << body.m << ","
            << "\"x\":" << body.x << ","
            << "\"y\":" << body.y << ","
            << "\"z\":" << body.z << ","
            << "\"vx\":" << body.vx << ","
            << "\"vy\":" << body.vy << ","
            << "\"vz\":" << body.vz
            << "}";
        if (i != bodies.size() - 1) oss << ",";
    }
    oss << "]";
    return oss.str();
}

Simulation simulation_;


/// OLAS
#define IP_ADDRESS "127.0.0.1"

class Server
{
public:
  struct CMDConnData
  {
    bool login;
  };
  using WSServer = websocket::WSServer<Server, CMDConnData>;
  using WSConn = WSServer::Connection;

  void run() {
    if (!wsserver.init(IP_ADDRESS, 1234)) {
      std::cout << "wsserver init failed: " << wsserver.getLastError() << std::endl;
      return;
    }

    running = true;
    ws_thr = std::thread([this]() {
      while (running.load(std::memory_order_relaxed)) {
        wsserver.poll(this);
        std::this_thread::yield();
      }
    });

    std::cout << "Server running..." << std::endl;
    ws_thr.join();
    std::cout << "Server stopped..." << std::endl;
  }

  void stop() { running = false; }


  // called when a new websocket connection is about to open
  // optional: origin, protocol, extensions will be nullptr if not exist in the request headers
  // optional: fill resp_protocol[resp_protocol_size] to add protocol to response headers
  // optional: fill resp_extensions[resp_extensions_size] to add extensions to response headers
  // return true if accept this new connection
  bool onWSConnect(WSConn& conn, const char* request_uri, const char* host, const char* origin, const char* protocol,
                   const char* extensions, char* resp_protocol, uint32_t resp_protocol_size, char* resp_extensions,
                   uint32_t resp_extensions_size) {
    struct sockaddr_in addr;
    conn.getPeername(addr);
    std::cout << "ws connection from: " << inet_ntoa(addr.sin_addr) << ":" << ntohs(addr.sin_port) << std::endl;
    std::cout << "request_uri: " << request_uri << std::endl;
    std::cout << "host: " << host << std::endl;
    if (origin) {
      std::cout << "origin: " << origin << std::endl;
    }
    if (protocol) {
      std::cout << "protocol: " << protocol << std::endl;
    }
    if (extensions) {
      std::cout << "extensions: " << extensions << std::endl;
    }
    return true;
  }

  // called when a websocket connection is closed
  // status_code 1005 means no status code in the close msg
  // status_code 1006 means not a clean close(tcp connection closed without a close msg)
  void onWSClose(WSConn& conn, uint16_t status_code, const char* reason) {
    std::cout << "ws close, status_code: " << status_code << ", reason: " << reason << std::endl;
  }

  // onWSMsg is used if RecvSegment == false(by default), called when a whole msg is received
  void onWSMsg(WSConn& conn, uint8_t opcode, const uint8_t* payload, uint32_t pl_len) {
    //std::cout << "onWSMsg: \n";

    if (opcode == websocket::OPCODE_PING) {
      conn.send(websocket::OPCODE_PONG, payload, pl_len);
      return;
    }
    if (opcode != websocket::OPCODE_TEXT) {
      conn.close(1003, "not text msg");
      return;
    }

// Bullshit text parse
    const char* data = (const char*)payload;
    const char* data_end = data + pl_len;
    char buf[4096] = {0};
    const char* argv[4096];
    char* out = buf + 1;
    int argc = 0;
    bool in_quote = false;
    bool single_quote = false;
    while (data < data_end) {
      char ch = *data++;
      if (!in_quote) {
        if (ch == ' ') *out++ = 0;
        else {
          if (*(out - 1) == 0) argv[argc++] = out;
          if (ch == '\'')
            in_quote = single_quote = true;
          else if (ch == '"')
            in_quote = true;
          else if (ch == '\\')
            *out++ = *data++;
          else
            *out++ = ch;
        }
      }
      else {
        if (single_quote) {
          if (ch == '\'')
            in_quote = single_quote = false;
          else
            *out++ = ch;
        }
        else {
          if (ch == '"')
            in_quote = false;
          else if (ch == '\\' && (*data == '\\' || *data == '"'))
            *out++ = *data++;
          else
            *out++ = ch;
        }
      }
    }
    if (argc) {
      *out = 0;
      std::string resp = onCMD(conn.user_data, argc, argv);
      if (resp.size()) conn.send(websocket::OPCODE_TEXT, (const uint8_t*)resp.data(), resp.size());
    }
  }

  // onWSSegment is used if RecvSegment == true, called when a segment is received
  // pl_start_idx: index in the whole msg for the 1st byte of payload
  // fin: whether it's the last segment
  void onWSSegment(WSConn& conn, uint8_t opcode, const uint8_t* payload, uint32_t pl_len, uint32_t pl_start_idx,
                   bool fin) {
    std::cout << "error: onWSSegment should not be called" << std::endl;
  }

 void broadcast(const std::string& msg) {
    wsserver.broadcast(msg);
  }


private:
  std::string onCMD(CMDConnData& conn, int argc, const char** argv) {

    std::string resp;
    std::string message(reinterpret_cast<const char*>((char *)argv[0]), argc);
    if (!strcmp(argv[0],"add")) {
        // Generate positions between -5.0 to +4.9 AU
        double x_au = (rand() % 100) / 10.0 - 5.0; 
        double y_au = (rand() % 100) / 10.0 - 5.0;

        // Convert AU to meters
        double x = x_au * 1.496e11;
        double y = y_au * 1.496e11;

        // Calculate distance from the Sun (assuming Sun is at origin)
        double r = std::sqrt(x*x + y*y);

        // Gravitational parameters (SI units)
        const double G = 6.67430e-11;      // Gravitational constant
        const double M_sun = 1.989e30;     // Sun's mass

        // Calculate required orbital velocity (circular orbit)
        double v_mag = std::sqrt(G * M_sun / r);

        // N-body simulation
        v_mag=10.0;
        x = x_au*1000000000.0;
        y = y_au*1000000000.0;
        // Remove to here to get sun orbit


        // Optional: Add slight randomness to velocity magnitude (0.9-1.1x)
        v_mag *= 0.9 + (rand() % 2001)/10000.0;  // Random factor between 0.9-1.1

        // Set velocity direction (tangential, counter-clockwise)
        double vx = -v_mag * y / r;
        double vy = v_mag * x / r;

        // N-body
        vx=10.0*(-10 + (rand() % 20));
        vy=10.0 *(-10 + (rand() % 20));

        std::cout << "Adding body at (" << x << ", " << y << ") with velocity (" << vx << ", " << vy << ")" << std::endl;


        // Create and add body (mass remains as before)
        double earth_mass = 0.5 * 1e24;
        Body body{earth_mass, x, y,0.0, vx, vy, 0.0, 0.0, 0.0};
        simulation_.addBody(body);
    }
    else if (!strcmp(argv[0], "remove")) {
        // Rmove body furthest from the origin
        simulation_.removeMostDistant();
        // simulation_.removeBody(0);
    }
    else if (!strcmp(argv[0], "force")) {
        simulation_.forceBody(0, 0.0, 0.0, 0.0);
    } else {
      resp = "invalid cmd";
    }
    return resp;

  }

private:
  WSServer wsserver;
  std::thread ws_thr;
  std::atomic<bool> running;
};


Server server;


int main() {    
    std::queue<std::string> msgQueue;
    std::mutex mtx;
    std::atomic<bool> running{true};


    // Start a simulation/broadcast thread
    std::thread simThread([&](){
        while (running) {
            // Step the simulation
            for (int i = 0; i < 30000; ++i) {
                simulation_.step(4.0);
            }

            // Serialize
            auto bodies = simulation_.getBodies();
            std::string json = serializeBodies(bodies);

            // Broadcast
            server.broadcast(json);

            // Sleep a bit so we don’t spam
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    });
    #if 0
    exec::static_thread_pool pool(4);
    // stdexec::scheduler auto sched = stdexec::thread_pool().scheduler();
    auto sched = pool.get_scheduler();

    std::thread execThread([&](){
        const double G = 6.67430e-11;  // Gravitational constant
        while (running) {
            // Step the simulation
            for (int i = 0; i < 30000; ++i) {

                // Reset accelerations
                for (auto& body : bodies_) {
                    body.ax = body.ay = body.az = 0.0;
                }
                double dt = 4.0;

                // Calculate forces in parallel for each i
                // 
                std::vector<stdexec::sender auto> force_senders;
                for (size_t i = 0; i < bodies_.size(); ++i) {                    
                    force_senders.emplace_back(
                        stdexec::just(i) | 
                        stdexec::then([this, G](size_t i) {
                            for (size_t j = i + 1; j < bodies_.size(); ++j) {
                                Body& a = bodies_[i];
                                Body& b = bodies_[j];
                                double dx = b.x - a.x;
                                double dy = b.y - a.y;
                                double dz = b.z - a.z;
                                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                                double F = G * a.m * b.m / (dist*dist + 1e-10);
                                double Fx = F * dx / dist;
                                double Fy = F * dy / dist;
                                double Fz = F * dz / dist;

                                a.ax += Fx / a.m;
                                a.ay += Fy / a.m;
                                a.az += Fz / a.m;
                                b.ax -= Fx / b.m;
                                b.ay -= Fy / b.m;
                                b.az -= Fz / b.m;
                            }
                        }) | stdexec::on(sched)
                    );
                }

                stdexec::sync_wait(stdexec::when_all(std::move(force_senders))).value();

                // Update velocities and positions
                for (auto& body : bodies_) {
                    body.vx += body.ax * dt;
                    body.vy += body.ay * dt;
                    body.vz += body.az * dt;
                    body.x += body.vx * dt + 0.5 * body.ax * dt * dt;
                    body.y += body.vy * dt + 0.5 * body.ay * dt * dt;
                    body.z += body.vz * dt + 0.5 * body.az * dt * dt;
                }
            }

            // Serialize
            auto bodies = simulation_.getBodies();
            std::string json = serializeBodies(bodies);

            // Broadcast
            server.broadcast(json);

            // Sleep a bit so we don’t spam
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    });
    #endif





    server.run();


#if 0
    // Deepseek suggestion
    // Using stdexec's thread pool scheduler
    exec::static_thread_pool pool(3);
     // Get a handle to the thread pool:
    auto sched = pool.get_scheduler();


    auto ws_thread = stdexec::schedule(sched) | stdexec::then([&] {
        while (running) {
#if 0
            server.poll(&handler);
            std::string msg;
            {
                std::lock_guard<std::mutex> lock(mtx);
                if (!msgQueue.empty()) {
                    msg = msgQueue.front();
                    msgQueue.pop();
                }
            }
            if (!msg.empty()) handler.broadcast(msg);
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
#endif
        }
    });

    auto sim_thread = stdexec::schedule(sched) | stdexec::then([&] {
        while (running) {
            sim.step(0.01);
            auto bodies = sim.getBodies();
            std::string json = serializeBodies(bodies);
            {
                std::lock_guard<std::mutex> lock(mtx);
                msgQueue.push(json);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    });

    // Wait for user input to exit
    std::cout << "Press Enter to exit..." << std::endl;
    std::cin.get();
    running = false;

    stdexec::sync_wait(stdexec::when_all(ws_thread, sim_thread));

#endif
    return 0;
}
