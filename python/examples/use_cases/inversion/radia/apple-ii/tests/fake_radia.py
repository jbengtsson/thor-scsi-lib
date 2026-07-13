from __future__ import annotations


class FakeRadia:
    def __init__(self):
        self.next_id = 1
        self.blocks = {}
        self.containers = {}
        self.transforms = []
        self.draw_attributes = []

    def _id(self):
        value = self.next_id
        self.next_id += 1
        return value

    def ObjRecMag(self, center, dimensions, magnetization):
        object_id = self._id()
        self.blocks[object_id] = {
            "center": list(center),
            "dimensions": list(dimensions),
            "magnetization": list(magnetization),
        }
        return object_id

    def ObjDrwAtr(self, object_id, color, line_width):
        self.draw_attributes.append((object_id, list(color), line_width))
        return object_id

    def ObjCnt(self, object_ids):
        object_id = self._id()
        self.containers[object_id] = list(object_ids)
        return object_id

    def TrfTrsl(self, vector):
        return ("translation", tuple(vector))

    def TrfOrnt(self, object_id, transform):
        self.transforms.append((object_id, transform))
        return object_id
