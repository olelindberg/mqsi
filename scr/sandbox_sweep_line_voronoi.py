import heapq
import math
import matplotlib.pyplot as plt

class Point:
    def __init__(self, x, y): self.x, self.y = x, y
    def __lt__(self, other): return self.y < other.y or (self.y == other.y and self.x < other.x)

class Edge:
    def __init__(self, start, direction): self.start, self.direction, self.end = start, direction, None

class Event:
    def __init__(self, point, arc=None): self.point, self.arc, self.valid = point, arc, True
    def __lt__(self, other): return self.point.y < other.point.y or (self.point.y == other.point.y and self.point.x < other.point.x)

class Arc:
    def __init__(self, site, prev=None, next=None):
        self.site, self.prev, self.next, self.event, self.edge = site, prev, next, None, None

class Voronoi:
    def __init__(self, sites):
        self.sites = sorted(sites, reverse=True)  # Max-heap
        self.event_queue = []
        for s in self.sites:
            heapq.heappush(self.event_queue, Event(s))
        self.arcs = None
        self.edges = []
        self.ly = 0

    def process(self):
        while self.event_queue:
            event = heapq.heappop(self.event_queue)
            self.ly = event.point.y
            if event.arc:
                if event.valid:
                    self.remove_arc(event)
            else:
                self.insert_arc(event.point)

    def insert_arc(self, p):
        if not self.arcs:
            self.arcs = Arc(p)
            return

        i = self.arcs
        while i:
            intersect, _ = self.intersect(p, i)
            if intersect:
                if i.next:
                    i.next.prev = Arc(i.site, prev=i, next=i.next)
                    i.next = i.next.prev
                else:
                    i.next = Arc(i.site, prev=i)
                i.next = Arc(p, prev=i, next=i.next)
                i.next.prev.next = i.next
                self.check_circle(i)
                self.check_circle(i.next)
                return
            i = i.next

        i = self.arcs
        while i.next:
            i = i.next
        i.next = Arc(p, prev=i)

    def remove_arc(self, event):
        arc = event.arc
        if arc.prev: arc.prev.next = arc.next
        if arc.next: arc.next.prev = arc.prev
        if arc.edge: arc.edge.end = event.point
        if arc.prev and arc.prev.edge: arc.prev.edge.end = event.point
        if arc.next and arc.next.edge: arc.next.edge.end = event.point
        self.check_circle(arc.prev)
        self.check_circle(arc.next)

    def check_circle(self, arc):
        if not arc or not arc.prev or not arc.next: return
        a, b, c = arc.prev.site, arc.site, arc.next.site
        if (b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y) >= 0:
            return
        A = b.x - a.x
        B = b.y - a.y
        C = c.x - a.x
        D = c.y - a.y
        E = A * (a.x + b.x) + B * (a.y + b.y)
        F = C * (a.x + c.x) + D * (a.y + c.y)
        G = 2 * (A * (c.y - b.y) - B * (c.x - b.x))
        if G == 0: return
        cx = (D * E - B * F) / G
        cy = (A * F - C * E) / G
        o = Point(cx, cy)
        r = math.hypot(a.x - cx, a.y - cy)
        p = Point(cx, cy - r)
        arc.event = Event(p, arc)
        heapq.heappush(self.event_queue, arc.event)

    def intersect(self, p, arc):
        if arc.site.y == p.y: return False, None
        a = 1 / (2 * (arc.site.y - self.ly))
        b = -2 * arc.site.x / (2 * (arc.site.y - self.ly))
        c = (arc.site.x**2 + arc.site.y**2 - self.ly**2) / (2 * (arc.site.y - self.ly))
        x = p.x
        y = a * x**2 + b * x + c
        return p.y < y, Point(x, y)

    def plot(self):
        for e in self.edges:
            if e.end:
                plt.plot([e.start.x, e.end.x], [e.start.y, e.end.y], 'k-')
        xs = [s.x for s in self.sites]
        ys = [s.y for s in self.sites]
        plt.plot(xs, ys, 'ro')
        plt.gca().set_aspect('equal')
        plt.show()

# Example usage:
sites = [Point(200, 200), Point(400, 200), Point(300, 400), Point(300, 100)]
v = Voronoi(sites)
v.process()
v.plot()
